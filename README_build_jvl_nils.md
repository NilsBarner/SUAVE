# README — Build & Rebuild JVL (MSYS2 MINGW64)

Concise, complete step-by-step guide for editing and rebuilding **JVL2.16** on Windows using the MSYS2 MINGW64 environment.

> Tested workflow (what you have): MSYS2 MINGW64 terminal, `gfortran`, `mingw32-make`, OpenBLAS installed.

---

## Quick summary

1. Install required MSYS2 packages (one-time).
2. Extract JVL source (avoid Drela's broken symlink).
3. Obtain/build `plotlib` (Xplot) separately.
4. Build `eispack`, then JVL.
5. To recompile after edits: edit sources in `src/` and run `mingw32-make` from `bin/` (Make will recompile only changed files).

---

## 0) Prerequisites (one-time)

Open a **MINGW64** shell (not MSYS):

```bash
# Update package DB + core system (may require restarting the shell)
pacman -Syu

# Install compilers, make, OpenBLAS and Ghostscript (for PS->PDF)
# Adjust packages to taste; these are the essentials used in this guide
pacman -S --needed mingw-w64-x86_64-toolchain mingw-w64-x86_64-make mingw-w64-x86_64-openblas mingw-w64-x86_64-ghostscript
```

---

## 1) Extract source safely

The tarball contains macOS symlinks. Extract everything **except** the problematic `plotlib` symlink:

```bash
cd /c/Users/nmb48/Downloads
# extract most files but exclude the bundled symlinked plotlib
tar --no-same-owner --no-same-permissions --warning=no-unknown-keyword -xvf jvl2.16.tgz --exclude='JVL2.16/plotlib'
# change to the top-level source
cd JVL2.16
ls
```

You should now see `bin/`, `src/`, `eispack/`, `runs/`, `bugs/`, etc.

---

## 2) Obtain `plotlib` (Xplot) and place in `JVL2.16/plotlib`

`plotlib` (aka Xplot plotting library) is distributed separately. Two routes:

**A. Preferred (download Xplot / copy from AVL source):**

* Download a copy of Drela's Xplot or copy the `plotlib` directory from a compatible AVL/XFOIL source tree.
* Put that directory at `JVL2.16/plotlib` (it must be a *real directory*, not a symlink).

**B. If you have the tarball with symlink and want a brute-force approach:**

* Extract the whole archive with `--dereference` on a Unix host or obtain a separate `plotlib` tree as above.

After you have a real `plotlib` directory, list it:

```bash
cd /c/Users/nmb48/Downloads/JVL2.16/plotlib
ls
# expected: Makefile, src/, ps/, lib/ etc.
```

---

## 3) Build `eispack` (linear algebra support)

```bash
cd /c/Users/nmb48/Downloads/JVL2.16/eispack
# use the mingw Makefile for Windows
mingw32-make -f Makefile.mingw
```

This produces `libeispack.a` in the `eispack/` directory.

---

## 4) Build `plotlib` (required for plotting)

`plotlib` supports multiple targets. JVL is built in **double precision** on MinGW. From `plotlib/`:

```bash
cd /c/Users/nmb48/Downloads/JVL2.16/plotlib
# Use the MinGW double-precision target and tell sub-makes to use mingw32-make
mingw32-make mingwDP MAKE=mingw32-make
```

Successful result: `libPlt_gDP.a` appears in `plotlib/`.

If the Makefile tries to call `make` recursively (no `make` exists), the `MAKE=mingw32-make` override fixes that.

---

## 5) Build JVL (final link step)

From the `bin/` directory run the provided Makefile:

```bash
cd /c/Users/nmb48/Downloads/JVL2.16/bin
mingw32-make -f Makefile.mingw
```

If everything is in place, the final compilation/link line will include:

```
../plotlib/libPlt_gDP.a ../eispack/libeispack.a -lopenblas
```

and `jvl.exe` will be produced in `JVL2.16/bin`.

Verify:

```bash
ls -l jvl.exe
./jvl.exe     # prints banner / prompt
```

---

## 6) Edit source & recompile (the usual workflow)

1. Open the source files you want to modify in `JVL2.16/src/` (e.g. `jvl.f`, `jmake.f`, `jaero.f`).
2. Save changes.
3. Re-run the Make from the `bin/` directory — **Make will recompile only changed files**:

```bash
cd /c/Users/nmb48/Downloads/JVL2.16/bin
# recommended: parallel build (replace 4 with your CPU cores)
mingw32-make -f Makefile.mingw -j4
```

### Useful variants

* Clean and rebuild from scratch (useful when module signatures change):

```bash
mingw32-make -f Makefile.mingw clean
mingw32-make -f Makefile.mingw
```

* If you changed `plotlib` sources, rebuild it too:

```bash
cd ../plotlib
mingw32-make mingwDP MAKE=mingw32-make
cd ../bin
mingw32-make -f Makefile.mingw
```

* If you modified `eispack` rebuild it before linking JVL:

```bash
cd ../eispack
mingw32-make -f Makefile.mingw
cd ../bin
mingw32-make -f Makefile.mingw
```

---

## 7) Debugging & common problems

**Missing `libPlt_gDP.a` at link:**
→ You must build `plotlib` (see step 4).

**`cannot find -lopenblas` at link:**
→ Install OpenBLAS (`pacman -S mingw-w64-x86_64-openblas`). Ensure you built JVL with the same MINGW64 environment.

**Recursive make fails with `make` not found:**
→ pass `MAKE=mingw32-make` to the top-level `mingw32-make` when building plotlib.

**Extraction fails due to Drela's symlink:**
→ extract archive excluding `JVL2.16/plotlib` then add a real `plotlib` directory.

**If you see many Fortran warnings about old constructs:**
→ these are expected (EISPACK and plotlib use legacy Fortran). Warnings are normal; only errors stop the build.

**If link fails with `ld` errors about missing symbols:**
→ paste the last ~20 lines of the linker output (they contain the real cause). Typically either a missing library or an incompatible object.

---

## 8) Tips for comfortable development

* Use version control:

```bash
cd /c/Users/nmb48/Downloads/JVL2.16
git init
git add -A
git commit -m "JVL2.16 base"
```

Commit small, incremental changes to be able to revert.

* Keep a small `build_jvl.sh` helper in the project root:

```bash
#!/usr/bin/env bash
set -e
cd $(dirname "$0")
# build libs
cd eispack && mingw32-make -f Makefile.mingw || exit 1
cd ../plotlib && mingw32-make mingwDP MAKE=mingw32-make || exit 1
# build jvl
cd ../bin && mingw32-make -f Makefile.mingw -j4
```

Make it executable with `chmod +x build_jvl.sh` and run `./build_jvl.sh`.

* For debugging, set lower optimization and keep symbols. Edit `Makefile.mingw` and change compiler flags (e.g. `-O2` → `-g -O0`).

---

## 9) How to produce and view plots

* JVL (with plotlib) produces PostScript files (`*.ps`). Convert to PDF with Ghostscript:

```bash
# assumed output.ps created by JVL
ps2pdf output.ps output.pdf
# or using ghostscript directly
gs -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile=output.pdf output.ps
```

Open the PDF with your preferred viewer.

---

## 10) If things go wrong — what to paste when asking for help

Always paste:

1. The **exact command** you ran.
2. The **last 20–40 lines** of terminal output (esp. linker errors).
3. `ls -l` output of the directory where you expected the library/binary (e.g. `JVL2.16/bin`).

That lets people diagnose link-time issues quickly.

---

If you want, I can also create a ready-to-use `build_jvl.sh` script (tailored to your current folders) and a minimal `debug` Make target. Want me to add those files to the README?
