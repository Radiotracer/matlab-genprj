mex   -DWIN32 '-IC:\mip\include' '-LC:\mip\lib64' -llibmiputil.lib -llibcl.lib -llibIRL.lib ... 
      -llibfftw3-3.lib -llibfftw3f-3.lib -llibfft-fftw3.lib -llibim.lib -llibimgio.lib  ...
     genprj.c GenprjSetup.c GetImages.c MeasToModPrj.c

 
 

 
clear;close all;
timg=readim('truth.im');
%figure;imshow(timg,[]);
prj=genprj('genprj.par',timg);
size(prj)
figure;imshow(squeeze(prj),[]);title('matlab genprj'); 
truth=readim('prj.im');
size(truth)
figure;imshow(squeeze(truth),[]);title('truth');
dif=truth-prj;
figure;imshow(squeeze(dif),[]);title('diff'); 





 imgprj=readim('prj.im');
 figure;imshow(squeeze(imgprj),[]);title('readim');
 size(imgprj)

t=imgprj-prj;
figure;imshow(squeeze(t),[]);
max
  
  
  mex   -DWIN32 '-IC:\mip\include' '-LC:\mip\lib64' -llibmiputil.lib -llibcl.lib -llibIRL.lib ... 
      -llibfftw3-3.lib -llibfftw3f-3.lib -llibfft-fftw3.lib -llibim.lib -llibimgio.lib  ...
     osem.c setup.c GetImages.c MeasToModPrj.c saveitercheck.c
 
 
 
prj=readim('prj.im');
figure;imshow(squeeze(prj),[]);title('readim');
size(imgprj)
recon=osem('osem.par',prj);
