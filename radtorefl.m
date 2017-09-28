% code to read in modtran-pcrtm clear sky and all sky results,
% in channel radiance form,
% solar zenith angle, and irradiances
% then convert to reflectance
%
% Input:
%     irad.mat - static irradiances and frequencies
%     list_tape5*.asc - lists multiple input tape5 filenames used as input to the PCRTM code,
%                      includes an index of the number of files contained
%     radiance*.bin - channel radiances produced by the PCRTM code for all tape5 runs
%
% Output:
%     reflectance_*.nc - netCDF file containing solar zenith angle and reflectances per PCRTM run

%     Read static irradiances and frequencies (nm)
%     Used to convert radiance to reflectance
      load('irad.mat');
      irad = new_irad;
      mx=new_x;
      mxchan = size(mx,1);

%     Determine how many list files were used and radiance files produced
%     (how many runs of the PCRTM)
      filelist=dir('list_tape5*.asc');
      radlist=dir('radiance*.bin');
      nrun=length(radlist);   % Number of times PCRTM was called; number of radiance files
      for ifovs=1:nrun;
%        Create netcdf output file
         fout = ['reflectance_' int2str(ifovs) '.nc']
         if (exist(fout))
            status = unix(['rm ' fout]);
         end
         nccreate(fout,'frequency','Dimensions',{'NumberChannels',mxchan},'Datatype','single')
         ncwrite(fout,'frequency',mx)

         idrun = fopen(filelist(ifovs).name,'r');
         idrad = fopen(radlist(ifovs).name,'r');
%        Read the number of FOVs in the PCRTM run from the list file
         nlist = fscanf(idrun,'%d');  % Number of FOVs run per one PCRTM run; or per list file
         
         save_sz = zeros(nlist,1);
         save_rf = zeros(mxchan,nlist);
         for ilist = 1:nlist;
            nchAll = fread(idrad,1,'double',24);
            solzen = fread(idrad,1,'double',24);
            rad = fread(idrad, nchAll,'double',24);

%           Calculate reflectance   
            sec=solzen./180*pi; 
            sec=cos(sec);
            fac=pi./(irad*sec);
            mref=rad.*fac;
            save_sz(ilist) = solzen;
            save_rf(1:nchAll,ilist) = mref;
         end
         nccreate(fout,'solzen','Dimensions',{'NumberFOVs',nlist})
         ncwrite(fout,'solzen',save_sz)
         nccreate(fout,'reflectance','Dimensions',{'NumberChannels',mxchan,'NumberFOVs',nlist},'Datatype','single')
         ncwrite(fout,'reflectance',save_rf)

%        ncdisp(fout)
      end
