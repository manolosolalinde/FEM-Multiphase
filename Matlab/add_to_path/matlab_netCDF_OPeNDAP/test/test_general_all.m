% test_general_all is called by test_all. According to the value of the
% variable test_type it will try to read a netcdf file or an opendap file. It
% calls the routine run_one_test many times with different values of the
% variable test_no determining exactly what reading test will be carried
% out.

% $Id: startup-xemacs.el,v 1.2 2003/10/17 04:46:07 man133 Exp man133 $
% Copyright J. V. Mansbridge, CSIRO, Fri Jun 16 17:38:37 EST 2006

% Use the test of ddsnc as a check of gross problems.

test_no = 0;
try
  run_one_test
catch
  disp(lasterr)
  return
end

% Test for a lot of other things that could go wrong.

num_tests = 117;
failure_list = [];
t_st_general = now;
for test_no = 1:num_tests
  % print out a message every 10 seconds to let the user know that the test
  % has not hung.
  
  t_now = now;
  elapsed_time = (t_now - t_st_general)*86400;
  if elapsed_time > 10
    disp([num2str(100*((test_no - 1)/num_tests)) '% done at' ...
	  datestr(t_now, 14)])
    t_st_general = t_now;
  end
  
  % Build up the failure list
  
  try
    run_one_test
  catch
    disp(lasterr)
    failure_list = [failure_list test_no];
  end
end

% Print out a message about how the tests went.

disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
num_failures = length(failure_list);
if num_failures == 0
  disp(['Passed all ' num2str(num_tests) ' tests'])
else
  disp(['There were ' num2str(num_failures) ' failures from ' ...
	num2str(num_tests) ' tests'])
  % Check whether all of the failures were when the variables were arrays of
  % characters.
  
  if all(ismember(failure_list, index_list_char))
    disp('All of the failures occurred reading variables that are arrays of')
    disp('characters. This is known to happen with the opendap implementation')
    disp('in Windows boxes. Since not many people use this type of variable')
    disp('it may not cause too much trouble')
  else
    disp('Something is seriously wrong - test_no for the failing tests was:')
    disp([num2str(failure_list)])
  end
end
