% SuperVoxel multi-modal segmentation

% Introduction
% The algorithm used in this code is the implementation of our paper [1]
% which was fulfilled in the Lab of Vision Engineering at the University of
% Lincoln (UK). It is a modification of the method Simple Linear Iterative
% Clustering (SLIC) which was proposed by Achanta et al. [2].
% Our method is optimized for medical images such as MRI, CT, etc. The
% contributions of our codes compared to conventional 2D and 3D superpixel
% are as follows:
% •	Multi-modal input (works for single-modal, as well)
% •	Taking the spatial resolution of the medical images into account, i.e.
% the voxel resolution in X and Y directions and the slice thickness.

clc
clear
close all

%% Initialization

% Load the parametrs ofrm th file "Opts.m"
run('Opts.m');
Integer = floor(Compactness);
Fraction = Compactness-Integer;
Cmpt = [num2str(Integer),'_',num2str(Fraction)];

% Data and results addresses
Case = '01';
ProtocolList = {'FLAIR','T1_Contrast','T2_Weighted','DTI_P_map','DTI_Q_map'}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InputImage_Path = [cd,'\Data\SGUL_Normalized'];
Output_Path = [cd,'\Results\SLIC'];
mkdir(Output_Path)






Image3D_All = [];
for P = 1:numel(ProtocolList)
    
    % Selecc the protocol
    Protocol = ProtocolList{P};
    
    % Displaying the patient case and protocols
    disp(['Case:  ',Case,'  Protocol: ',Protocol])
    
    % Load the data
    ImageFileName = ['SGUL_Normalized_Case_',Case,'_',Protocol,'.mat'];
    load(fullfile(InputImage_Path,ImageFileName));
    Image3D = Image3D_Normalized;
    Image3D_temp = Image3D;
    Min_I = min(Image3D_temp(:));
    Max_I = max(Image3D_temp(:));
    Image3D_temp = (Image3D_temp-Min_I)/(Max_I-Min_I);
    eval(['I_',num2str(P),' = Image3D_temp;']);
end

% Initial arrays
I = I_1;
S = nthroot((voxel_X * voxel_Y * voxel_Z),3);
X_window = floor(size(I,1)/voxel_X);
Y_window = floor(size(I,2)/voxel_Y);
Z_window = floor(size(I,3)/voxel_Z);
r1 = ceil(voxel_X/2);
r2 = ceil(voxel_Y/2);
r3 = ceil(voxel_Z/2);
[X,Y,Z] = meshgrid((2:Y_window-1).*voxel_X-r2,(2:X_window-1).*voxel_Y-r1,(2:Z_window-1).*voxel_Z-r3);

% Initial Centers for 5 protocols
for P = 1:numel(ProtocolList)
    eval(['I_C_',num2str(P),' = zeros(size(X));']);
end

for i=1:size(X,1)
    for j = 1:size(X,2)
        for l = 1:size(X,3)
            eval( ['I_C_',num2str(P),' (i,j,l) = I_',num2str(P),'(X(i,1,1), Y(1,j,1),Y(1,1,l));'] );
        end
    end
end

C = [X(:),Y(:),Z(:)];
for P = 1:numel(ProtocolList)
    eval(['C = cat(2,C,I_C_',num2str(P),'(:));']);
end

k = size(C,1);

% Update supervoxel centers
for i = 1:k
    temp_C = C(i,:);
    x = temp_C(1);
    y = temp_C(2);
    z = temp_C(3);
    
    temp_W = I_1(x-1:x+1,y-1:y+1,z);
    a = abs(gradient(temp_W));
    [x,y] = find(a==min(a(:)), 1, 'first');
    C(i,1) = temp_C(1)+x-2;
    C(i,2) = temp_C(2)+y-2;
    
    for P = 1:numel(ProtocolList)
        eval(['C(i,3+',num2str(P),') = I_',num2str(P),'(C(i,1), C(i,2));']);
    end

end

% Assigning the lables
Label = -1*ones(size(I_1));
Distance = inf*ones(size(I_1));

counter = Iterations; % Iterastions
while counter
    for c_k = 1:k
        Center = [C(c_k,1),C(c_k,2),C(c_k,3)];
        for i = Center(1)-2*r1+1:Center(1)+2*r1-1
            for j = Center(2)-2*r2+1:Center(2)+2*r2-1
                for l = Center(3)-2*r3+1:Center(3)+2*r3-1
                    % Calculating the spatial distance
                    d_c_term = 0;
                    d_c_temp = 0;
                    for P = 1:numel(ProtocolList)
                        eval(['d_c_temp = (I_',num2str(P),'(Center(1),Center(2),Center(3))-I_',num2str(P),'(i,j,l))^2;']);
                        d_c_term = d_c_term + d_c_temp;
                    end
                    d_c = sqrt(d_c_term );
                    
                    % Calculating the intensity distance
                    d_s = sqrt((term_x*(Center(1)-i))^2 + (term_y*(Center(2)-j))^2+(term_z*(Center(3)-l))^2);
                    
                    % Calculating the total distance D between C_k and i
                    D = sqrt((d_c)^2+(d_s/S)^2*Compactness^2);
                    if D < Distance(i,j,l)
                        Distance(i,j,l) = D;
                        Label(i,j,l) = c_k;
                    end
                end % for l
            end % for j
        end % for i
    end % for SV Label
    
    %% Updating stage
    counter = counter-1;
    
    % Compute new cluster centers
    if counter > 0
        for c_k = 1:k
            [x1,y1,z1] = find3((Label==c_k));
            C(c_k,1) = round(mean(x1));
            C(c_k,2) = round(mean(y1));
            C(c_k,3) = round(mean(z1));
        end
        
        % make sure that the centers are inside image dimensions
        C(C(:,1)<2*r1,1) = 2*r1;
        C(C(:,1)>size(I,1)-2*r1,1) =size(I,1)-2*r1;
        C(C(:,2)<2*r2,2) = 2*r2;
        C(C(:,2)>size(I,2)-2*r2,2) =size(I,2)-2*r2;
        C(C(:,3)<2*r3,3) = 2*r3;
        C(C(:,3)>size(I,3)-2*r3,3) =size(I,3)-2*r3;
    end % update
end % While
SLIC_Labels_3D = Label;


%% Save
% Save the supervoxel map volumes into MAT file
Output_Name = fullfile(Output_Path,['MRI_SLIC_Labels_Size',num2str(voxel_X),...
    'x',num2str(voxel_Y),'x',num2str(voxel_Z),'_Compactness_0',Cmpt,'_Case_',num2str(Case),'.mat']);
save (Output_Name,'SLIC_Labels_3D');



