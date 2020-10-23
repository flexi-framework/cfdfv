# Fortran CI

In Github there is the possibility to run specific build commands every time a code change is pushed to the master branch. The file that controls what exact commands are run is
```
.github/workflows/fortran.yml
```
If you take look at that file, you will see in the first line the name of this 'Github Action'. Then follows the conditions for when the script is to be executed. The last block defines the jobs that are run and on which platform. You will see that there are three steps for this job

1. `make shared`, which builds the shared library
2. `make`, which builds the `cfdfv` executable
3. `make check`, which runs `cfdfc` on a specified set of test cases and compares the results with a given reference

# How are the checks performed?

In Addition to the `yml` file, there are two more shell scripts and one python file added in the `utils` directory.
```
utils/makeCheck.sh
utils/csvdiff.py
utils/updateReferences.sh
```
These files are used during the `make check` job, mentioned earlier. When you take a look at `makeCheck.sh` you will see that this script runs a specified test case, using the compiled binary `cfdfv` and compares the latest time step output to a reference file (using `cgnsdiff` and `csvdiff.py` when appropriate). The only thing that is checked is, if the values match the reference values (to a small tolerance). Because checking a test case also requires running that test case, not all test cases are active. To activate them, just remove the `#` in front of the test case that you wish to add.

If you wish to update the reference solutions, there is a separate make target called `make updateReferences`. Just execute this command from the root directory and all reference cases will be updated. The corresponding script that performs the updating (`updateReferences.sh`) also has most of the test cases commented out, so make sure that the cases that you want to update are active.
