# Setup: Running the project on the DCC machines (MPI)

Quick guide to connect to the machines, copy the project over, and run the benchmarks.

## 1. Connect via SSH (after following the professor's tutorial)

Replace `l815` with the machine you want.

```bash
ssh -i ~/.ssh/upXXXXYYYYY -J upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt upXXXXYYYYY@l815.alunos.dcc.fc.up.pt
```

The `-J` option uses `ssh.alunos.dcc.fc.up.pt` as a jump host to reach `l815`.

## 2. Set up SSH between machines (inside `l815`)

Generate a key **without a passphrase** (just press Enter when prompted):

```bash
ssh-keygen -t ed25519
```

Copy the key to another machine (probably not even necessary, since the home directory is shared):

```bash
ssh-copy-id l816
```

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
(otherwise `mpirun` may hang).
