clear classes;
reload = 1;
nhl = 1;
fields = 0;
plotCorrelations = 0;
includeKEmods = 1;
includeENmods = 1;
handFit = 0;
doFit = 1;
plotResults = 0;
useStart = 0;
pstart =  [-0.5 3 7 0 0 0];
params = 1:22; % bond lengths to include
envs = 1:100; % environments to include in fit

if (reload)
   load('h2bonds/h2BondDat.mat');
end

if (doFit || handFit)
   disp('building models');
   m = cell(1,size(params,2));
   for ipar = params
      m{ipar} = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});
   end
   if (includeKEmods)
      mixKEdiag = Mixer([0 0],2);
      mixKEbond = Mixer(0,1);
      for ipar = params
         m{ipar}.addKEmodDiag(1,1,mixKEdiag);
         m{ipar}.addKEmodBonded(1,1,1,1,mixKEbond);
      end
   end
   if (includeENmods)
      mixENdiag = Mixer([0 0],2);
      mixENbond = Mixer(0,1);
      for ipar = params
         m{ipar}.addENmodDiag(1,1,mixENdiag);
         m{ipar}.addENmodBonded(1,1,1,1,mixENbond);
      end
   end
end

if (doFit)
   disp('Starting to do parameter fitting');
   f1 = Fitme;
   for ipar = params
      f1.addFrag(m{ipar},HL{ipar});
   end
   f1.includeKE = includeKEmods;
   f1.includeEN = includeENmods * ones(1,6);
   f1.setEnvs(envs);
   nfitpar = f1.npar;
   start = zeros(1,nfitpar);
   if (useStart)
      start = pstart;
   end
   % code from using matlabs optimization toolbox
   %limits = 3 * ones(1,nfitpar);
   %options = optimset('DiffMinChange',1.0e-5);
   %pt = lsqnonlin(@f1.err, start,-limits,limits,options);
   options = LMFnlsq;
   options.Display =1;
   options.FunTol = 1.0e-6;
   options.XTol = 1.0e-5;
   [pfit, Ssq, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
end

if (plotResults)
   disp('Starting to do plots');
   if (doFit)
      pt = pfit;
   else
      pt = pstart;
   end
   ic = 0;
   for ipar = params
      ll = LL{ipar,1};
      hl = HL{ipar};
      t1 = ll.EEhf;
      rr = (ic+1):(ic+size(t1,2));
      lke(rr) = t1;
      hke(rr) = hl.EEhf;
      le1(rr) = ll.Een(1);
      he1(rr) = hl.Een(1);
      disp(['starting calc on ipar ',num2str(ipar)]);
      m{ipar}.setPars(pt);
      m{ipar}.solveHF;
      mke(rr) = m{ipar}.EEhf;
      me1(rr) = m{ipar}.Een(1);
      ic = ic + size(t1,2);
   end
   figure(100);
   hold off;
   plot(lke,hke,'r.');
   hold on;
   plot(lke,lke,'k.');
   plot(lke,mke,'b.');
   title('ke')
   figure(200);
   plot(mke,hke,'g.');
   title('ke')
   
   figure(101);
   hold off;
   plot(le1,he1,'r.');
   hold on;
   plot(le1,le1,'k.');
   plot(le1,me1,'b.');
   title('EN')
   figure(201);
   plot(me1,he1,'g.');
   title('EN')
end
