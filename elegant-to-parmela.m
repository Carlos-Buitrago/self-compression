clear all
mc2   = 510.99906E-3;					% e- rest mass [MeV]
%Ebar = mean(A(:,6))					% average e- energy in MeV
%gam  = Ebar/mc2;					% average Lorentz factor
fRF   = 2856E6;							% RF frequency [Hz]
c     = 2.99792458E8;					% speed of light [m/s]
lamRF = c/fRF;

filename = 'C:\Users\cfern\Documents\Tesi\Autocompressione\ASTRA_to_TStep\1pC-oncrest.001.sdds'
delimiterIn = ' ';
headerlinesIn = 10;
E = importdata(filename,delimiterIn,headerlinesIn);
E2T=E.data;
E2Tf(:,1)=1e2*E2T(:,1);
E2Tf(:,3)=1e2*E2T(:,3);
E2Tf(:,2)=E2T(:,2);
E2Tf(:,4)=E2T(:,4);
E2Tf(:,5)=-(360*fRF*(E2T(:,5)-mean(E2T(:,5))));  % multiplying everything by -1 to rotate
E2Tf(:,6)=(E2T(:,6)- mean(E2T(:,6)))*mc2;
N = length(E2T(:,6))
%CHANGE OUTPUT FILENAME HERE
filename2 = 'C:\Users\cfern\Documents\Tesi\Autocompressione\ASTRA_to_TStep\TStep_In'
fnout = [filename2 '\1pC-oncrest.txt'];
%if yn == 'y'												% write Tstep input txt file
  if exist(fnout)
    cmnd = ['delete ' fnout];
    eval(cmnd)
  end
  fid = fopen(fnout,'w');
  for j = 1:N
    str = [];
    str = [str sprintf('%11.9e ',E2Tf(j,1))];					% X [cm]
    str = [str sprintf('%11.9e ',E2Tf(j,2))];					% X' [rad]
    str = [str sprintf('%11.9e ',E2Tf(j,3))];					% Y [cm]
    str = [str sprintf('%11.9e ',E2Tf(j,4))];					% y' [rad]
    str = [str sprintf('%11.9e ',-E2Tf(j,5))];					% phi [deg]
    str = [str sprintf('%11.9e\n',E2Tf(j,6))];		% E/mc^2 [MeV]
    fprintf(fid,str);
  end
  disp(' ')
  disp(['Tstep input file written: ' fnout])
  fclose(fid);
%end
