# YALMIP + ECOS Setup (Windows, MATLAB R2025b, MinGW64)

Required for the YALMIP/ECOS variant in `HM2_powered_descent/main_task2.m`.
These instructions document the exact steps verified on this machine.

---

## 1. Download toolboxes (PowerShell)

```powershell
$tb = "$env:USERPROFILE\Documents\MATLAB\toolboxes"
New-Item -ItemType Directory -Force -Path $tb

# YALMIP
Invoke-WebRequest "https://github.com/yalmip/YALMIP/archive/refs/heads/master.zip" -OutFile "$tb\yalmip.zip"
Expand-Archive "$tb\yalmip.zip" -DestinationPath $tb -Force

# ecos-matlab wrapper
Invoke-WebRequest "https://github.com/embotech/ecos-matlab/archive/refs/heads/master.zip" -OutFile "$tb\ecos-matlab.zip"
Expand-Archive "$tb\ecos-matlab.zip" -DestinationPath $tb -Force

# ECOS C source (required for compilation — not included in the ZIP, git submodule)
git clone https://github.com/embotech/ecos.git "$tb\ecos-matlab-master\ecos"
```

## 2. Compile ECOS (MATLAB Command Window)

The standard `makemex` script fails on R2025b + MinGW64 for two reasons:
- The library is now called `libut.lib` (not `ut.lib`)
- MATLAB's linker incorrectly picks up the Windows Store Python stub path when resolving `-lut`

Workaround: compile components via `makemex`, then link manually.

```matlab
% Step 1 — compile object files (stops at linking, that is expected)
cd(fullfile(getenv('USERPROFILE'), 'Documents', 'MATLAB', 'toolboxes', 'ecos-matlab-master', 'bin'));
makemex   % will fail at "Linking..." — ignore the error

% Step 2 — link manually using the correct MinGW64 libut path
ut_lib = fullfile(matlabroot, 'extern', 'lib', 'win64', 'mingw64', 'libut.lib');
d = '-largeArrayDims -DDLONG -DLDL_LONG';
cmd = sprintf(['mex %s "%s" amd_1.obj amd_2.obj amd_aat.obj amd_control.obj ' ...
    'amd_defaults.obj amd_dump.obj amd_global.obj amd_info.obj amd_order.obj ' ...
    'amd_post_tree.obj amd_postorder.obj amd_preprocess.obj amd_valid.obj ' ...
    'ldl.obj kkt.obj preproc.obj spla.obj cone.obj ecos.obj ctrlc.obj ' ...
    'timer.obj splamm.obj equil.obj ecos_mex.obj ecos_bb.obj ' ...
    'ecos_bb_preproc.obj wright_omega.obj expcone.obj -output "ecos"'], d, ut_lib);
eval(cmd);
```

`MEX completed successfully` confirms the binary `ecos.mexw64` was created in the `bin/` folder.

## 3. Add paths to MATLAB startup (MATLAB Command Window)

Run once — persists across MATLAB restarts.

```matlab
startup_path = fullfile(getenv('USERPROFILE'), 'Documents', 'MATLAB', 'startup.m');
tb = fullfile(getenv('USERPROFILE'), 'Documents', 'MATLAB', 'toolboxes');
fid = fopen(startup_path, 'a');
fprintf(fid, "addpath(genpath('%s'));\n", fullfile(tb, 'YALMIP-master'));
fprintf(fid, "addpath('%s');\n",          fullfile(tb, 'ecos-matlab-master', 'bin'));
fclose(fid);
```

## 4. Verify

```matlab
% YALMIP version
yalmip('version')               % should print e.g. '20250626'

% ECOS binary (which takes precedence over ecos.m in the same folder)
which ecos                      % should point to ecos.mexw64

% End-to-end SOCP test
x = sdpvar(2,1);
res = optimize(norm(x) <= 1, x(1), sdpsettings('solver','ecos','verbose',0));
fprintf('Status: %d  x = [%.4f, %.4f]\n', res.problem, value(x(1)), value(x(2)));
% Expected: Status: 0  x = [-1.0000,  0.0000]
```

## Notes

- `exist('ecos','file')` returns `2` (finds `ecos.m` in the same folder) rather than `3`.
  This is a MATLAB path-ordering quirk; `which ecos` correctly resolves to the MEX binary,
  and YALMIP calls it correctly. The availability check in `main_task2.m` uses
  `exist(...) > 0` which evaluates to `true` for both 2 and 3.
- MinGW64 must be installed (MATLAB Add-Ons → MinGW-w64 C/C++ Compiler).
- Tested on: Windows 11, MATLAB R2025b, YALMIP 20250626, MinGW64 12.x.
