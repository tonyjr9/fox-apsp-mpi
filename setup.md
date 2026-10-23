# Setup: Running the project on the DCC machines (MPI)

How to connect to the Lab8 machines (`l801`–`l829`), copy the project
over, and run the benchmarks. Replace `upXXXXYYYYY` with your login
name.

## 1. Set up SSH access (once, on your own PC)

Generate an SSH key. The default file is `~/.ssh/id_ed25519`:

```sh
ssh-keygen -t ed25519 -C upXXXXYYYYY
```

Leave the passphrase **empty** — `mpirun` cannot type a passphrase and
hangs silently if the key has one. If `~/.ssh/id_ed25519` already
exists, keep it and skip the generation (overwriting it would break
everything else that uses it).

Copy the key to the jump host and to the machine in the students'
area. Both ask for your up password once:

```sh
ssh-copy-id upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt
ssh-copy-id -o ProxyJump=upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt \
  upXXXXYYYYY@kanopus.alunos.dcc.fc.up.pt
```

The home directory is shared across all Lab8 machines, so authorizing
on `kanopus` once covers `l801`–`l829`. No `-i` flag or
`~/.ssh/config` entries are needed afterwards: ssh offers the default
`~/.ssh/id_ed25519` on every hop.

## 2. Connect to a machine

```sh
ssh -J upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt upXXXXYYYYY@l815.alunos.dcc.fc.up.pt
```

`-J` uses the jump host to reach the machine (note the full name,
`l815.alunos.dcc.fc.up.pt` — it only resolves through the jump host).
Log in **at least once to every machine** that appears in your
hostfile so its host key is accepted; otherwise `mpirun` cannot reach
it. To logout run `exit`.

## 3. Check which machines are up

From inside a lab session:

```sh
for h in $(awk '{print $1}' hostfile.txt); do
  ssh -o BatchMode=yes -o ConnectTimeout=3 $h true 2> /dev/null \
    && echo "$h ok" || echo "$h DOWN"
done
```

A quick alternative is `ping $h` from inside a lab machine. Machines
go down for mundane reasons (a student switched the machine off).
Remove any DOWN machine from `hostfile.txt` before running, otherwise
`mpirun` may hang. If **every** machine shows DOWN, the SSH setup in
section 1 is broken, not the machines.

## 4. Copy the project to the server

Run this **on your own PC** (not on the server), from the folder that
contains `fox-apsp-mpi/`:

```sh
rsync -avz -e "ssh -J upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt" \
  fox-apsp-mpi/ upXXXXYYYYY@l815.alunos.dcc.fc.up.pt:~/Paralela/
```

## 5. The hostfile

One machine per line with `slots` = number of cores; `hostfile.txt` in
this repo is an example. **The first machine in the hostfile must be
the machine you are logged in to**, and write it exactly as
`hostname` prints it (the lab hosts report an uppercase name like
`L815`); a first entry that is not your machine, or that mismatches
in case, makes `mpirun` fail with `ORTE does not know how to route
a message`.

## 6. Run

### Silence the X11 authorization warning

On the lab machines, `mpirun` may print
`Authorization required, but no authorization protocol specified`.
This comes from hwloc (the hardware topology library used by Open MPI)
trying to open the X display of whoever is logged in at the console.
It is harmless, but disable that probe before running anything:

```sh
export HWLOC_COMPONENTS=-gl
```

It must be exported in the shell, so `mpirun` itself gets it. `unset
DISPLAY` or `mpirun -x HWLOC_COMPONENTS=-gl` do not help. To make it
permanent, add the line to the **top** of `~/.bashrc`, before any line
that returns early for non-interactive shells.

### Run

On the server, inside the project folder (the program takes the input
matrix file and an output file):

```sh
mpirun -n 32 --hostfile hostfile.txt fox input6.bin output6.bin
```

By default `mpirun` fills the first machine in the hostfile up to its
`slots` limit, and only then starts launching processes on the second
machine. To launch round robin (one process on each machine at a
time):

```sh
mpirun -np 32 --hostfile hostfile.txt --map-by node fox input6.bin output6.bin
```

### Run the benchmarks

```sh
NPS="1 4 16 32 64" MPIRUN="mpirun --hostfile hostfile.txt" ./bench/run_bench.sh
```

- `NPS` – the process counts to test.
- `MPIRUN` – the MPI command, here using `hostfile.txt` to distribute
  processes across machines.