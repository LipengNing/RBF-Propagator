function CollectIndices(input_dwi,mask_name,ExportName,flag)
if(~exist('flag','var'))
    flag=0;
end

if(~matlabpool('size'))
    matlabpool open;
end


if(flag==1)
% stronger constraints on the decay of signal, slower
[D,V,mask] = run_GaussianBasis_Decay(input_dwi,mask_name);
else
% weaker constraints, but faster, it may takes about 10 hours for one volume
[D,V,mask] = run_GaussianBasis(input_dwi,mask_name);
end

save([ExportName,'_coef']);

%% write data into Nift files 
nii=MRIread(mask_name);
 
RTOP=ComputeRTOP(D,V,mask);
RTOP=permute(RTOP,[2 1 3]);
nii.vol=RTOP;
nii.fspec=[ExportName,'_RTOP.nii'];
MRIwrite(nii,[ExportName,'_RTOP.nii']);
%%
[RTAP, RTPP]=ComputeRTAPandRTPP(D,V,mask);
RTAP=permute(RTAP,[2 1 3]);
RTPP=permute(RTPP,[2 1 3]);

nii.vol=RTAP;
nii.fspec=[ExportName,'_RTAP.nii'];
MRIwrite(nii,[ExportName,'_RTAP.nii']);

nii.vol=RTPP;
nii.fspec=[ExportName,'_RTPP.nii'];
MRIwrite(nii,[ExportName,'_RTPP.nii']);

%%
MSD=ComputeMSD(D,V,mask);
MSD=permute(MSD,[2 1 3]);
nii.vol=MSD;
nii.fspec=[ExportName,'_MSD.nii'];
MRIwrite(nii,[ExportName,'_MSD.nii']);
%%
MFD=ComputeMFD(D,V,mask);
MFD=permute(MFD,[2 1 3]);
nii.vol=MFD;
nii.fspec=[ExportName,'_MFD.nii'];
MRIwrite(nii,[ExportName,'_MFD.nii']);
%%
[HD, SDD]=ComputeHDandSDD(D,V,mask);
HD=permute(HD,[2 1 3]);
SDD=permute(SDD,[2 1 3]);
nii.vol=HD;
nii.fspec=[ExportName,'_HD.nii'];
MRIwrite(nii,[ExportName,'_HD.nii']);

nii.vol=SDD;
nii.fspec=[ExportName,'_SDD.nii'];
MRIwrite(nii,[ExportName,'_SDD.nii']);
%%
[GK,NIG]=ComputeKurtosis(D,V,mask);
GK=permute(GK,[2 1 3]);
NIG=permute(NIG,[2 1 3]);
nii.vol=GK;
nii.fspec=[ExportName,'_GK.nii'];
MRIwrite(nii,[ExportName,'_GK.nii']);

nii.vol=NIG;
nii.fspec=[ExportName,'_NIG.nii'];
MRIwrite(nii,[ExportName,'_NIG.nii']);

%%
QMSD=ComputeMSD_Qspace(D,V,mask);
QMSD=permute(QMSD,[2 1 3]);
nii.vol=QMSD;
nii.fspec=[ExportName,'_QMSD.nii'];
MRIwrite(nii,[ExportName,'_QMSD.nii']);
%

QMFD=ComputeMFD_Qspace(D,V,mask);
QMFD=permute(QMFD,[2 1 3]);

nii.vol=QMFD;
nii.fspec=[ExportName,'_QMFD.nii'];
MRIwrite(nii,[ExportName,'_QMFD.nii']);

%%
D=squeeze(D);
fa=zeros(1,size(D,3));
d_bar=zeros(1,size(D,3));
parfor i=1:size(D,3)
    [fa_i,d_i]=tensor2fa(D(:,:,i));
    fa(i)=fa_i;
    d_bar(i)=d_i;
end
Fa=zeros(size(mask));
Fa(mask~=0)=fa;
D_bar=zeros(size(mask));
D_bar(mask~=0)=d_bar;
Fa=permute(Fa,[2 1 3]);
D_bar=permute(D_bar,[2 1 3]);
nii.vol=Fa;
nii.fspec=[ExportName,'_Fa.nii'];
MRIwrite(nii,[ExportName,'_Fa.nii']);

nii.vol=D_bar;
nii.fspec=[ExportName,'_Trace.nii'];
MRIwrite(nii,[ExportName,'_Trace.nii']);

