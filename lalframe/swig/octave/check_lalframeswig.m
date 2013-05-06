## Check SWIG Octave module wrapping lalframe
## Author: Karl Wette, 2011, 2012

## check module load
lalframe;
assert(exist("lalframe", "var"));
assert(exist("lalframecvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
lalcvar.lalDebugLevel = bitor(LALERROR, LALMEMDBG);
disp("passed module load");

# check object parent tracking
a = lalframe.new_lalframeswig_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalframecvar.lalframeswig_test_parent_map.s;
  lalframecvar.lalframeswig_test_parent_map.s = lalcvar.lalswig_test_struct_const;
endfor
clear a b c;
CheckMemoryLeaks();
disp("passed object parent tracking");

## passed all tests!
disp("PASSED all tests");
