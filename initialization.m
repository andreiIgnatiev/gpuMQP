
global nc dt
%rr

%CHROMOSOME ADDITIONS
%%
PtsOnChr=10;%15 %# of points on chromsome blob structure
%%
n_chr=8;%46; %# of chromosomes
%%
L_chr=1;%0.75; %Total length of chromosome arms, used to create blob structure
%%
X_chr=3*rand(2,n_chr)-1.5; %Centers of chromosomes
KL = 0.30;
Krad = KL*ones(1,PtsOnChr);
delta = linspace(0,2*pi,PtsOnChr);
r=L_chr*ones(1,PtsOnChr)+0.5*rand(1,PtsOnChr)-0.25;
theta=linspace(0,2*pi-2*pi/(PtsOnChr),PtsOnChr);

alpha=zeros(1,n_chr);


for i=1:n_chr
    %Dimension of ChrBlobPts(x/y, 20 points on each blob, n_chr
    %chromosomes)
    KTPts(1,:,i)=2*Krad.*cos(delta)+X_chr(1,i);
    KTPts(2,:,i)=2*Krad.*sin(delta)+X_chr(2,i);
    %X-coordinates of blob of chromosomes
    ChrBlobPts(1,:,i)=2*r.*cos(theta)+X_chr(1,i);
    %Y-coordinates of blob of chromosomes
    ChrBlobPts(2,:,i)=0.5*r.*sin(theta)+X_chr(2,i);
    
end

% figure (26)
% plot(X_chr(1,:),X_chr(2,:),'or','MarkerSize',3,'LineWidth',1)
% hold on
% for i=1:n_chr
%     plot(ChrBlobPts(1,:,i),ChrBlobPts(2,:,i),'r','LineWidth',1)
% end



d_k= 0.34; %10 Mb in length of centromere %0.34
%how big the centromeres are lol 

r_k=(d_k/2)*ones(1,n_chr);



xRadius=15;      % changed 
yRadius=15;
X_loc=1*(xRadius/3)*randi([-1,1],1,nc).*rand(1,nc);
Y_loc=1*(yRadius/3)*randi([-1,1],1,nc).*rand(1,nc);

%N=20;
%vec=randperm(N);
%xLoc=vec(randperm(length(vec),1));
%yLoc=vec(randperm(length(vec),1));
%phi=rand(1)*(pi/2);
%X_loc=xLoc*cos(phi);
%Y_loc=yLoc*sin(phi);
%X_loc=[-2.5 2.5];
%Y_loc=[0 0];
%X_chr_L
%chr_centers = 
centers = [X_loc; Y_loc];
radii = 0.3*ones(1,nc);
%size of centrosomes
c_rad = radii(1); 



%number of MTs
n_MT = 300; %
%need to initialize a state
%centrosome attachments
c_attach = randi([1 nc],1,n_MT); %randi=creates random integer value between 1 and nc (#'s in brackets) and 1,n_MT is the size of the vector
cortex_S = zeros(1,n_MT);
cortex_B = zeros(1,n_MT);

ChrCentMTInteractions=zeros(2,n_MT);%2xn_MT, 1st row=chr #, 2nd row centrosome #. 0 if no interaction, # otherwise


%added chr things
chr_attachS = zeros(1,n_MT);
chr_attachB = zeros(1,n_MT);

%MT lengths
MT_length = 0.5*rand(1,n_MT);
%MT angle from centrosome
MT_angle = rand(1,n_MT)*2*pi; 
%Center of centrosome MT is attached to
for i = 1:n_MT
    mt_centers(:,i) = [centers(:,c_attach(i))];
end
%initial MT location
X_MT = [(MT_length+c_rad).*cos(MT_angle(1,:))+mt_centers(1,:); (MT_length+c_rad).*sin(MT_angle(1,:))+mt_centers(2,:)];
mt_vec(1,:)=(1./(MT_length)).*(X_MT(1,:)-(mt_centers(1,:)+c_rad*cos(MT_angle)));
            %unit vector -> scale or divide by lngth
            %x-location of mt minus x-location of peri centriolar material
mt_vec(2,:)=(1./(MT_length)).*(X_MT(2,:)-(mt_centers(2,:)+c_rad*sin(MT_angle)));

%Elliptical boundary of cell
xCenter = 0;
yCenter = 0;
%xRadius = 15;
%yRadius = 15;
%domain for plotting
xmin = -xRadius-1;
xmax = xRadius+1;
ymin = -yRadius-1;
ymax = yRadius+1;
zeta = linspace(0,2*pi,1000);
x = xRadius * cos(zeta) + xCenter;
y = yRadius * sin(zeta) + yCenter;
%plot(x, y, 'LineWidth', 3);
%%axis equal;
%xlim([-xRadius-1 xRadius+1]);
%ylim([-yRadius-1 yRadius+1]);
%hold on

%States
n_state = 4;
MT_state = ones(1,n_MT);
checkstate1=size(find(MT_state==1));
checkstate2=size(find(MT_state==2));
checkstate3=size(find(MT_state==3));
checkstate4=size(find(MT_state==4));
%state 1: growing; %state 2: shrinking; %state 3: slipping (on cortex);
%state 4: binding (to cortex); state 5: slipping on chromosome arms;
%state 6: binding to chromosome arms; state 7: bound to kinetochore
