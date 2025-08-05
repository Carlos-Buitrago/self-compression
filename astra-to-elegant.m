clear all

%%%costanti
mc2   = 510.99906E-3;					% e- rest mass [MeV]
fRF   = 2856E6;							% RF frequency [Hz]
c     = 2.99792458E8;				% speed of light [m/s]
cathode=0; %%per .ini
% comb=1;
sdds_need=1;
TStep=0;

lamRF = c/fRF;
% q=6e-9;


dzslice=10e-6;

need=0;
% q=75e-12;
slice_anal=1;
long_generation=0;

% CHANGE INPUT FILENAME HERE
ffn = 'C:\Users\cfern\Documents\Tesi\Autocompressione\ASTRA_to_TStep'
       ffn=[ffn '\1pC-oncrest.001'];

if exist(ffn)

% ffn
delimiterIn = ' ';
% headerlinesIn =1;%540075;%3060411;%
E = importdata(ffn,delimiterIn);
%(E(2:end,3)+E(1,3),E(2:end,2)*1e3,'o')
end

    [a,b]=hist(E(2:end,3)/c,98);
    A_AS(:,5)=-E(2:end,3)-mean(E(2:end,3));   %%z[m]
A_AS(:,11)=E(2:end,3)/c-mean(E(2:end,3)/c); %t[s]




%%%unit conversion to be checked %%%

%% E_ASTRA_sq eV^2
% E_ASTRA_sq=((mc2*1e6)^2+(E(2:end,4).^2+E(2:end,5).^2+(E(2:end,6)+E(1,6)).^2));
% pippo=1e6*mc2+(-2+sqrt(4+4*(E(2:end,4).^2+E(2:end,5).^2+(E(2:end,6)+E(1,6)).^2)))/2;
%
% % E_ASTRA=sqrt(E_ASTRA_sq);
% [a,b]=hist(1e-6*((E_ASTRA)-mean(E_ASTRA)),180);

p_mod_sq=1e-6^2*(E(2:end,4).^2+E(2:end,5).^2+(E(2:end,6)+E(1,6)).^2);
% E_s=sqrt(mc2^2+p_mod_sq)+mc2;

E_s=sqrt(mc2^2+p_mod_sq);
[h,s]=hist(E_s-mean(E_s),15);

if long_generation==1
subplot(121)
plot(b,a,'r')
subplot(122)
plot(s,h)
end
% [a,b]=hist(1e-6*((pippo)-mean(pippo)),178);
%
% plot(b,a,'r')

%%%transverse

if exist(ffn)
A_AS(:,[1 3])=E(2:end,[1 2]); %%x[m] y[m]
% A_AS(:,[2 4])=E(:,[2 4]);
A_AS(:,2)=E(2:end,4)./(sqrt(E(2:end,4).^2+E(2:end,5).^2+(E(2:end,6)+E(1,6)).^2)); %%x' [rad]
A_AS(:,4)=E(2:end,5)./(sqrt(E(2:end,4).^2+E(2:end,5).^2+(E(2:end,6)+E(1,6)).^2)); %%y' [rad]
% A_AS(:,5)=E(2:end,3)-mean(E(2:end,3));   %%z[m]
A_AS(:,6)=(E_s-mean(E_s))/mean(E_s);   %%dE/E
A_AS(:,8)=E_s;                             %E[MeV]
% A_AS(:,11)=E(2:end,3)/c-mean(E(2:end,3)/c); %t[s]
A_AS(:,12)=E(2:end,8);
gam=mean(E_s)/mc2;
E_bar=mean(E_s);
S      = cov(A_AS);				    % 6X6 beam covariance matrix
ex     = sqrt(det(S(1:2,1:2)));		% rms x-emittance [m-rad]
ey     = sqrt(det(S(3:4,3:4)));	    % rms y-emittance [m-rad]
gex    = gam*ex;					% rms norm x-emittance [m-rad]
gey    = gam*ey;					% rms norm y-emittance [m-rad]
sigx   = sqrt(S(1,1));				% sigx [m]
sigy   = sqrt(S(3,3));				% sigy [m]
betax  = S(1,1)/ex;					% betax [m]
alphax = -S(1,2)/ex;				% alphax [m]
betay  = S(3,3)/ey;					% betay [m]
alphay = -S(3,4)/ey;				% alphax [m]
sz     = std(A_AS(:,5));				% rms bunch length [m]
sd     = std(A_AS(:,6));				% rms rel. E-spread [unitless]






end

if slice_anal==1
         A=A_AS;


if need==0
       prova1=A;
end
prova1=sortrows(prova1,5);


dz0=(max(prova1(:,5))-min(prova1(:,5)));
zslice_min=min(prova1(:,5));
zslice_max=zslice_min+dzslice;
n=round(dz0/dzslice);
nslice=1;
%n=1;


    for nslice=1:n
        i=1;
        k=1;
        prova2_exist=0;
        clear prova2


    for i=1:length(prova1(:,5))

        if (prova1(i,5)<zslice_max && prova1(i,5)>zslice_min )
           prova2(k,:)=prova1(i,:);
           k=k+1;
           prova2_exist=1;
        end

    end


        if prova2_exist==1 && length(prova2(:,1))>2
           Ebar_nslice=mean(prova2(:,8));
           prova2(:,10) = (prova2(:,8) - Ebar_nslice)/Ebar_nslice;
           S      = cov(prova2);				    % 6X6 beam covariance matrix
           gam_1(nslice,1)  = Ebar_nslice/mc2;
           ex_1(nslice,1)     = sqrt(det(S(1:2,1:2)));		% rms x-emittance [m-rad]
           ey_1(nslice,1)     = sqrt(det(S(3:4,3:4)));	    % rms y-emittance [m-rad]
           gex_1(nslice,1)    = gam_1(nslice,1)*ex_1(nslice,1);					% rms norm x-emittance [m-rad]
           gey_1(nslice,1)     = gam_1(nslice,1)*ey_1(nslice,1);					% rms norm y-emittance [m-rad]
           betax_1(nslice,1)    = S(1,1)/ex_1(nslice,1)  ;					% betax [m]
           alphax_1(nslice,1)   = -S(1,2)/ex_1(nslice,1)  ;				% alphax [m]
           betay_1(nslice,1)    = S(3,3)/ey_1(nslice,1)  ;					% betay [m]
           alphay_1(nslice,1)   = -S(3,4)/ey_1(nslice,1)  ;				% alphax [m]
           st_1(nslice,1)     = std(prova2(:,11));				% rms bunch length [m]
           sd_1 (nslice,1)    = std(prova2(:,10));				% rms rel. E-spread [unitless]
           q_slice            = -sum(prova2(:,12))*1e-9;
           I_slice(nslice,1)=(c*q_slice)/(dzslice);
%          I_slice(nslice,1)=(c*q/length(A))*length(prova2(:,1))/(dzslice);

           zn_slice(nslice,1)=mean(prova2(:,5))*1e6;
           zslice_mean(nslice,1)=zslice_min+dzslice;

        end

        zslice_min=zslice_max;
        zslice_max=zslice_min+dzslice;





    end

ind = find(sum(gex_1,2)==0) ;
gex_1(ind,:) = [] ;
ind = find(sum(gey_1,2)==0) ;
gey_1(ind,:) = [] ;
ind = find(sum(sd_1,2)==0) ;
sd_1(ind,:) = [] ;
ind = find(sum(I_slice,2)==0) ;
I_slice(ind,:) = [] ;
ind = find(sum(zslice_mean,2)==0) ;
zslice_mean(ind,:) = [] ;
ind = find(sum(betax_1,2)==0) ;
betax_1(ind,:) = [] ;
ind = find(sum(alphax_1,2)==0) ;
alphax_1(ind,:) = [] ;



figure(9)

subplot(131)
plot(zslice_mean(:,1)*1E6,(gex_1(:,1))*1e6,'b')
hold on
% plot(zslice_mean(:,1)*1E6,gey_1(:,1)*1e6,'r')
% plot(zslice_mean(:,1)*1E6,sd_1(:,1)*1E2,'m')
% plot(zslice_mean(:,1)*1E6,I_slice(:,1)*1e-2,'k')

ylabel('{\epsilon_n} [mm-mrad]','FontSize',20)
xlabel('{\itz} [\mum]','FontSize',20)
subplot(132)

% hold on
plot(zslice_mean(:,1)*1E6,betax_1(:,1),'r')
%  plot(zslice_mean(:,1)*1E6,alphax_1(:,1),'--b')
hold on
plot(zslice_mean(:,1)*1E6,alphax_1(:,1),'--r')
plot(zslice_mean(:,1)*1E6,(I_slice(:,1)),'k')
hold on
ylabel('{I_{slice}} [A]','FontSize',20)
xlabel('{\itz} [\mum]','FontSize',20)
%legend('{\DeltaE/<E>} [%]',' {I_{FWHM}} [kA]')
%
%
  subplot(133)
% plot(zslice_mean(:,1)*1E6,betax_1(:,1),'b')
  plot(zslice_mean(:,1)*1E6,sd_1(:,1)*1E2,'m')
hold on
% plot(zslice_mean(:,1)*1E6,alphax_1(:,1),'r')
% ylabel('{I_{slice}} [A]','FontSize',20)
xlabel('{\itz} [\mum]','FontSize',20)
legend('\deltaE/E [%]','\alpha_x')
% subplot(122)

%  figure(7)
%
% subplot(121)
%
% plot(zslice_mean(:,1)*1E6,betax_1(:,1),'b')
% subplot(122)
% ylabel('{\beta_x} [m]','FontSize',20)
% xlabel('{\itz} [\mum]','FontSize',20)
% % hold on
% % plot(zslice_mean(:,1)*1E6,betay_1(:,1),'r')
% % plot(zslice_mean(:,1)*1E6,alphax_1(:,1),'--b')
% % hold on
% % plot(zslice_mean(:,1)*1E6,alphay_1(:,1),'--r')
% plot(zslice_mean(:,1)*1E6,sd_1(:,1)*1E2,'m')
% hold on
% ylabel('{energy spread} [%]','FontSize',20)
% xlabel('{\itz} [\mum]','FontSize',20)
% %legend('{\DeltaE/<E>} [%]',' {I_{FWHM}} [kA]')
% %
%
% subplot(122)
% hist(prova2(:,5)*1e6)
% hold on
% hist(A(:,5)*1e6)
% % plot(zslice_mean(nslice,1)*1E6,gex_1(nslice,1),'.b')
% % hold on
% % plot(zslice_mean(nslice,1)*1E6,gey_1(nslice,1),'.r')
% % ylabel('{\epsilon_n} [mm-mrad]','FontSize',20)
% % xlabel('{\itz} [\mum]','FontSize',20)
%
%
%
end
%%%%%%plot%%%%%%%%

if slice_anal==0
   A=A_AS;
end

figure(10)
%hold off
subplot(221)
plot(A(:,1)*1E3,A(:,2)*1E3,'.b')
hold on
%plot_ellipse(inv(1E6*9*S(1:2,1:2)),1,'r-')
xlabel('{\itx} [mm]','FontSize',20)
ylabel('{\itx}\prime [mrad]','FontSize',20)
title(['{\it\gamma\epsilon_x}=' sprintf('%6.2f',1E6*gex) ' \mum'],'FontSize',20)
%enhance_plot('times',14,1,1)
hold on
subplot(222)
plot(A(:,3)*1E3,A(:,4)*1E3,'.g')
hold on
%plot_ellipse(inv(1E6*9*S(3:4,3:4)),1,'r-')
xlabel('{\ity} [mm]','FontSize',20)
ylabel('{\ity}\prime [mrad]','FontSize',20)
title(['{\it\gamma\epsilon_y}=' sprintf('%6.2f',1E6*gey) ' \mum'],'FontSize',20)
%enhance_plot('times',14,1,1)
hold on

subplot(223)
plot(A(:,5)*1E3,A(:,6)*1E2,'.r')
hold on
xlabel('{\itz} [mm]','FontSize',20)
ylabel('\Delta{\itE}/\langle{\itE}\rangle [%]','FontSize',20)
title(['{\it\sigma_z}=' sprintf('%3.0f',1E6*sz) ' \mum, {\it\sigma_E}/\langle{\itE}\rangle=' sprintf('%6.2f',1E2*sd) '%'],'FontSize',20)

N=i;
v = axis;
h = v(4) - v(3);
text(2*v(2),v(4)-0*h/7,['{\itN}=' sprintf('%8.0f',N)])
text(2*v(2),v(4)-1*h/7,['\langle{\itE}\rangle=' sprintf('%7.3f',E_bar) ' MeV'])
text(2*v(2),v(4)-2*h/7,['{\it\sigma_x}=' sprintf('%7.3f',sigx*1E3) ' mm'])
text(2*v(2),v(4)-3*h/7,['{\it\beta_x}=' sprintf('%7.3f',betax) ' m'])
text(2*v(2),v(4)-4*h/7,['{\it\alpha_x}=' sprintf('%7.3f',alphax)])
text(2*v(2),v(4)-5*h/7,['{\it\sigma_y}=' sprintf('%7.3f',sigy*1E3) ' mm'])
text(2*v(2),v(4)-6*h/7,['{\it\beta_y}=' sprintf('%7.3f',betay) ' m'])
text(2*v(2),v(4)-7*h/7,['{\it\alpha_y}=' sprintf('%7.3f',alphay)])
%enhance_plot('times',14,1,1)




if sdds_need==1

fnout = [ffn '.sdds'];
%if yn == 'y'												% write Elegant input ASCII file
  if exist(fnout)
    cmnd = ['delete ' fnout];
    eval(cmnd)
  end
  fid = fopen(fnout,'w');
  fprintf(fid,'SDDS1\n');
  fprintf(fid,'&column name=x, units=m, type=double,  &end\n');
  fprintf(fid,'&column name=xp, symbol=x'', type=double,  &end\n');
  fprintf(fid,'&column name=y, units=m, type=double,  &end\n');
  fprintf(fid,'&column name=yp, symbol=y'', type=double,  &end\n');
  fprintf(fid,'&column name=t, units=s, type=double,  &end\n');
  fprintf(fid,'&column name=p, units="m$be$nc", type=double,  &end\n');
  fprintf(fid,'&data mode=ascii, &end\n');
  fprintf(fid,'! page number 1\n');
  fprintf(fid,'                %5.0f \n',N);
  for j = 1:N
    str = [];
    str = [str sprintf('%23.21e ',A(j,1))];					% X [m]
    str = [str sprintf('%23.21e ',A(j,2))];					% X' [rad]
    str = [str sprintf('%23.21e ',A(j,3))];					% Y [m]
    str = [str sprintf('%23.21e ',A(j,4))];					% y' [rad]
    str = [str sprintf('%23.21e ',A(j,5)/c)];					% Z [m] (NOT ROTATED)
    str = [str sprintf('%23.21e\n',(A(j,6)+1)*E_bar/mc2)];		% E/mc^2 [gamma]
    fprintf(fid,str);
  end
  disp(' ')
  disp(['Elegant input file written: ' fnout])
  fclose(fid);
end

% % % index_mid = round(N_step/2);
% % % beta_x0 = beta_x(index_mid);
% % % alfa_x0 = alfa_x(index_mid);
% % % mismatch_param = 1/2*(beta_x0*gamma_x-2*alfa_x0*alfa_x+gamma_x0*beta_x);

% plot(A_AS(:,1)*1e3,A_AS(:,2)*1e3,'.r')
% plot(A_AS(:,5)*c*1e3,A_AS(:,6)*1e2,'.')
%
% [a1,b1]=hist(E(2:end,4)./(E(2:end,6)+E(1,6)));
% hold on
% plot(b1,a1)
% [a2,b2]=hist(TS(:,2));
% plot(b2,a2)
