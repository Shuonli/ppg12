# Building yaml-cpp for PPG12

The ROOT macros read their YAML configs through yaml-cpp
(`#include <yaml-cpp/yaml.h>` in 59 tracked `.C` files). The sPHENIX `new`
release does not ship it (`$OFFLINE_MAIN/lib*` and `$OPT_SPHENIX/lib*` have no
`libyaml-cpp*`), so it lives in the private prefix `$MYINSTALL`
(`/sphenix/u/shuhang98/install`), next to the anatreemaker library
`lib/libCaloAna24.so`. This page records exactly how the installed library was
built and how to rebuild it into another prefix.

## Installed library (as of 2026-09-10)

| item | value |
|---|---|
| source | `https://github.com/jbeder/yaml-cpp.git`, commit `73ef0060aaa1a9dc742c1d8a36fa336b35e94035` (`git describe`: `0.8.0-70-g73ef006`, 2024-12-30, "Avoid including \<iostream\> in library code") |
| local checkout | `/sphenix/user/shuhangli/yaml-cpp/yaml-cpp` at that commit, with no modified tracked files (the only `git status` entry is the untracked `build_shared/` directory, `build/` is gitignored). Build tree `build/` (2025-04-04) is the one that produced the installed files. `build_shared/` (2025-01-17, gcc 12.1.0, SL7) is an older build of the same source and is **not** what is installed. |
| toolchain | `source /opt/sphenix/core/bin/sphenix_setup.sh -n new` on Alma 9: g++ 14.2.0 (`/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/gcc/14.2.0-2f0a0/x86_64-el9/bin/g++`), cmake 3.31.0, Unix Makefiles |
| cmake cache | `CMAKE_INSTALL_PREFIX=/sphenix/u/shuhang98/install`, `CMAKE_INSTALL_LIBDIR=lib64`, `BUILD_SHARED_LIBS=ON` (command line), `YAML_BUILD_SHARED_LIBS=ON`, `YAML_CPP_BUILD_TESTS=OFF`, `YAML_CPP_BUILD_TOOLS=ON`, `YAML_CPP_BUILD_CONTRIB=ON`, `CMAKE_BUILD_TYPE` empty |
| compile flags | `-std=gnu++11 -fPIC -Wall -Wextra -Wshadow -Weffc++ -Wno-long-long -pedantic -pedantic-errors` (yaml-cpp's own CMakeLists pins C++11. The ABI that matters is libstdc++ from the same gcc 14.2.0 that built ROOT 6.32.06 in the `new` release, requiring `GLIBCXX_3.4.29` / `CXXABI_1.3`) |
| output | `lib64/libyaml-cpp.so.0.8.0` (1577992 bytes, md5 `6a5a0aa3e27337696183743a3569e121`, SONAME `libyaml-cpp.so.0.8`) plus symlinks `libyaml-cpp.so.0.8` and `libyaml-cpp.so`, `include/yaml-cpp/`, `lib64/cmake/yaml-cpp/*.cmake`, `lib64/pkgconfig/yaml-cpp.pc` (Version 0.8.0). 45 files in `build/install_manifest.txt`. |
| leftover | `lib64/libyaml-cpp.a` (2025-04-04 15:29, nine minutes before the `.so`) is in no install manifest and was compiled with GCC 12.1.0 (the SL7 toolchain, `strings -a libyaml-cpp.a \| grep GCC:`), so it did not come from the gcc 14.2.0 `build/` tree. Its origin is not recorded and nothing in PPG12 links it. |

## Recipe

Run from any directory with write access. `MYINSTALL` is the target prefix.
Use the literal `/sphenix/u/shuhang98/install` layout (`lib64/`) if you need the
existing macros to work unchanged, see the caveat below.

```bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n new      # gcc 14.2.0, cmake 3.31.0
export MYINSTALL=/sphenix/u/shuhang98/install              # or your own prefix

git clone https://github.com/jbeder/yaml-cpp.git
git -C yaml-cpp checkout 73ef0060aaa1a9dc742c1d8a36fa336b35e94035

cmake -S yaml-cpp -B yaml-cpp/build \
      -DCMAKE_INSTALL_PREFIX="$MYINSTALL" \
      -DCMAKE_INSTALL_LIBDIR=lib64 \
      -DBUILD_SHARED_LIBS=ON \
      -DYAML_CPP_BUILD_TESTS=OFF
cmake --build   yaml-cpp/build -j8
cmake --install yaml-cpp/build
```

Tag `0.8.0` (`f7320141120f720aecc4c32be25586e7da9eb978`) is 70 commits older than
the pinned commit. It builds and links the same way but is not byte-identical.

Rebuilding with this recipe from the local checkout into a scratch prefix
(2026-09-10) reproduced `libyaml-cpp.so.0.8.0` byte for byte (same md5, same
4008 exported symbols). Byte identity depends on the source path: with
`CMAKE_BUILD_TYPE` empty there is no `-DNDEBUG`, so the `assert()` messages embed
the absolute source path (12 strings in the installed `.so`). A fresh clone at
the same commit in another directory gives a functionally identical library
(same 4008 exported symbols) with a different md5.

## Verify

```bash
source /sphenix/user/shuhangli/ppg12/env.sh          # puts $MYINSTALL/lib64 on LD_LIBRARY_PATH
                                                     # and $MYINSTALL/include on ROOT_INCLUDE_PATH
root -l -b -q -e 'gSystem->Load("libyaml-cpp.so"); std::cout << gSystem->GetLibraries() << std::endl;'
root -l -b -q -e '#include <yaml-cpp/yaml.h>' \
              -e 'gSystem->Load("libyaml-cpp.so"); std::cout << YAML::LoadFile("/sphenix/user/shuhangli/ppg12/efficiencytool/config_bdt_nom.yaml")["output"]["var_type"].as<std::string>() << std::endl;'
```

Expected: the library path printed by `GetLibraries()` is under `$MYINSTALL/lib64`
(for the default prefix ROOT prints the automount form
`/direct/sphenix+u/shuhang98/install/lib64/libyaml-cpp.so.0.8.0`) and the second
command prints `bdt_nom`. The `#include` has to be its own `-e` argument, a
multi-line `-e` does not parse.

## Caveat: the macros still hardcode the path

59 tracked macros, including the whole nominal chain
(`efficiencytool/RecoEffCalculator_TTreeReader.C:47`, `MergeSim.C:8`,
`merge_periods.C:147`, `CalculatePhotonYield.C:84`,
`plotting/plot_final_selection.C:17`, `FunWithxgboost/BDTinput.C:31`), call

```cpp
gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
```

with the literal path, and three tracked docs repeat it as house convention
(`.claude/rules/root-macros.md`, `.claude/agents/code-writer.md`,
`wiki/guides/common-issues.md`). Until those lines are migrated to
`gSystem->Load("libyaml-cpp.so")` (resolved through the `LD_LIBRARY_PATH` that
`env.sh` sets, tested 2026-09-10) or to `$MYINSTALL`, the literal path must exist
with `lib64/libyaml-cpp.so` inside it, whatever `MYINSTALL` points at.
