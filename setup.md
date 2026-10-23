# Setup: Running the project on the DCC machines (MPI)

Quick guide to connect to the machines, copy the project over, and run the benchmarks.

## 1. Connect via SSH (after following the professor's tutorial)

Replace `l815` with the machine you want.

```bash
ssh -i ~/.ssh/upXXXXYYYYY -J upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt upXXXXYYYYY@l815.alunos.dcc.fc.up.pt
```

The `-J` option uses `ssh.alunos.dcc.fc.up.pt` as a jump host to reach `l815`.

## 2. Set up SSH between machines (inside `l815`)

The lab machines only accept key login (no passwords), and `mpirun` can't type a passphrase,
so the key used between machines **must not have a passphrase**. If your existing key has one,
`mpirun` hangs silently waiting for it. Create a separate key with an empty passphrase:

```bash
ssh-keygen -t ed25519 -N "" -f ~/.ssh/id_ed25519_lab
```

Authorize it. The home directory is shared, so doing this once covers every machine:

```bash
cat ~/.ssh/id_ed25519_lab.pub >> ~/.ssh/authorized_keys
chmod 700 ~/.ssh
chmod 600 ~/.ssh/authorized_keys
```

SSH ignores `authorized_keys` if `~/.ssh` or the file can be written by others, hence the `chmod`.

Tell ssh to use this key for the lab machines by adding to `~/.ssh/config`:

```
Host l8*
    IdentityFile ~/.ssh/id_ed25519_lab
```

Check that it works without any prompt (it should exit silently):

```bash
ssh -o BatchMode=yes l816 true
```

If it prints `Permission denied (publickey)`, ssh is still offering a key with a passphrase or a
key not listed in `authorized_keys`.

### Accept the host keys of all machines

MPI needs to connect to each machine without asking for confirmation, so you have to connect
**at least once** to each one. Instead of using the long command above, from inside the SSH session you can simply run:

```bash
ssh l801 hostname
```

Repeat this for every machine in `hostfile.txt` (`l801`, `l802`, ...).

## 3. Copy the project to the server

Run this **on your own PC** (not on the server), from the folder that contains `fox-apsp-mpi/`:

```bash
rsync -avz -e "ssh -J upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt" fox-apsp-mpi/ upXXXXYYYYY@l815.alunos.dcc.fc.up.pt:~/Paralela/
```

Alternative: `git clone` on the server, but this requires adding the SSH key to GitHub (annoying).

## 4. Run the benchmarks

### Silence the X11 authorization warning

On the lab machines, `mpirun` may print
`Authorization required, but no authorization protocol specified`. This comes from hwloc (the
hardware topology library used by Open MPI) trying to open the X display of whoever is logged in
at the console. It is harmless, but disable that probe before running anything:

```bash
export HWLOC_COMPONENTS=-gl
```

It must be exported in the shell, so `mpirun` itself gets it. `unset DISPLAY` or
`mpirun -x HWLOC_COMPONENTS=-gl` do not help. To make it permanent, add the line to the
**top** of `~/.bashrc`, before any line that returns early for non-interactive shells.

### Run

On the server, inside the project folder:

```bash
NPS="1 4 16 32 64" MPIRUN="mpirun --hostfile hostfile.txt" ./bench/run_bench.sh
```

- `NPS` – the process counts to test.
- `MPIRUN` – the MPI command, here using `hostfile.txt` to distribute processes across machines.

## 5. Check which machines are up or down

```bash
for h in $(awk '{print $1}' hostfile.txt); do
  ssh -o BatchMode=yes -o ConnectTimeout=3 $h true 2>/dev/null && echo "$h ok" || echo "$h DOWN"
done
```

If any machine shows `DOWN`, remove it from `hostfile.txt` before running the benchmark
(otherwise `mpirun` may hang). If **every** machine shows `DOWN`, the problem is almost certainly
the SSH key (see section 2), not the machines.
