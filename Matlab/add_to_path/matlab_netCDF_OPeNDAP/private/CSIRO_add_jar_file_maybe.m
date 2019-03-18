function CSIRO_add_jar_file_maybe

% CSIRO_add_jar_file_maybe might add toolsUI.jar to the dynamic java path.

option = 1;
switch option
 case 1
  % Assume that toolsUI.jar is already in the path (static or
  % dynamic) and so do nothing.
  
  return
 case 2
  % Try to add toolsUI.jar to the dynamic path. Note that the use
  % of addpath will clear result in clear all being run and any
  % global variables will be lost.
  
  % If there is a problem with the java virtual machine (like matlab
  % running with the -nojvm option) then don't try to change the
  % dynamic java path.

  if ~isempty(javachk('jvm'))
    return
  end

  % Predict the full name of the jar file assuming that it is in the same
  % directory as getnc.

  directory = which('getnc');
  path_name = directory(1:(end-7));
  name_jar_file = [path_name 'toolsUI.jar'];

  % If the jar file exists then check whether it is already in the
  % java path. If it is not there then add it.

  if exist(name_jar_file)
    jpath = javaclasspath('-all');
    not_there = 1;
    for ii = 1:length(jpath)
      if strcmp(name_jar_file, jpath{ii})
	not_there = 0;
	break
      end
    end
    if not_there
      javaaddpath(name_jar_file)
    end
  end
end
