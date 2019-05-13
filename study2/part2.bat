@echo off

set v0=1
set tw=1.5708
set dT=0.01
set stageOrTime=1
set log_stages=0


REM ----------
REM r0 = 1
REM ----------

..\sim.exe 0.01 %v0% 1 %tw% 12000000 %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 0.25 %v0% 1 %tw% 1000000  %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 0.5  %v0% 1 %tw% 400000   %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 1    %v0% 1 %tw% 200000   %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 1.5  %v0% 1 %tw% 140000   %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 2    %v0% 1 %tw% 100000   %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 2.5  %v0% 1 %tw% 80000    %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 3    %v0% 1 %tw% 60000    %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 3.5  %v0% 1 %tw% 60000    %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 4    %v0% 1 %tw% 40000    %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 4.5  %v0% 1 %tw% 40000    %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
..\sim.exe 5    %v0% 1 %tw% 40000    %dT% a e 1 0.7 %stageOrTime% %log_stages% 0
