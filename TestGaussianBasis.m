function YorN=TestGaussianBasis()

load TestData;


[D,V] = GaussianBasis_Coef(S,u,b);
mask=ones(1,9);
YorN=1;

%% Test RTOP
RTOP_true=[2792        3268        3919        9473       10995        5924        5288        2632        3159];
RTOP=ComputeRTOP(D,V,mask);
if(sum(round(RTOP/1e2)~=RTOP_true))
    YorN=0;
    sprintf('RTOP is wrong!');
end
%% Test MSD
MSD_true=[372   350   215   133    26   173   166   376   354];
MSD=ComputeMSD(D,V,mask);
if(sum(round(MSD*1e6)~=MSD_true))
    YorN=0;
    sprintf('MSD is wrong!');
end
%% Test ODF
ODF_true=[    45    52    45
    53    55    56
    44    53    49
    56    74    71
     0     4     0
    20     4     1
    46    36    49
    58    52    54
    67    65    71
];

[ODF,gg]=ComputeODF(D,V,mask);
ODF=reshape(ODF,9,2562);ODF(ODF<0)=0;
if(sum(sum(round(ODF(:,1:3)*1000)~=ODF_true)))
    YorN=0;
    sprintf('ODF is wrong!');
end


