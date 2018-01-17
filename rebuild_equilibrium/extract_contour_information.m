clear contour_length contour_length_raw psi_value_list_raw X_contour_raw Z_contour_raw X_contour Z_contour X_contour_down Z_contour_down X_contour_up Z_contour_up
clear X_contour_shift Z_contour_shift psi_value_list length_contour_down length_contour_up

Ncontours=find_number_of_contours(Cont_psi2D)
if Ncontours~=Nradial
    disp('WARNING: the number of flux surfaces does not match the number of contours !!!!');
end
psi_value_list_raw=zeros(Ncontours,1);
contour_length_raw=zeros(Ncontours,1);

% X1_Nradial(1)=round(mid_X+(X_axis/DX));
% X2_Nradial(1)=round(mid_X+(X_axis/DX));
index=0;

for (n=1:Ncontours)
    index=index+1;
    psi_value_list_raw(n)=Cont_psi2D(1,index);
    contour_length_raw(n)=Cont_psi2D(2,index);
    for s=1:contour_length_raw(n)
        index=index+1;
        X_contour_raw(n,s)=Cont_psi2D(1,index);
        Z_contour_raw(n,s)=Cont_psi2D(2,index);
    end
%     X1_Nradial(n)=round(mid_X+(min(X_contour(n,1:contour_length(n))))/DX);
%     X2_Nradial(n)=ceil(mid_X+(max(X_contour(n,1:contour_length(n))))/DX);
end

for (n=1:Ncontours)
    if contour_length_raw(n)<2
	    X_contour_raw(n,1:2)=[X_axis X_axis];
        Z_contour_raw(n,1:2)=[Z_axis Z_axis];
        contour_length_raw(n)=2;
	end
end

%[min_val SMALLEST_CONTOUR]=min(contour_length_raw(2:end))
%if SMALLEST_CONTOUR>1
%    X_contour_raw(SMALLEST_CONTOUR+1,1:2)=[X_axis X_axis];
%    Z_contour_raw(SMALLEST_CONTOUR+1,1:2)=[Z_axis Z_axis];
%    contour_length_raw(SMALLEST_CONTOUR+1)=2;
%    psi_value_list_raw(SMALLEST_CONTOUR+1)=psi_global;
%else
%    [min_val SMALLEST_CONTOUR]=min(contour_length_raw)
%    X_contour(1,1:2)=[X_axis X_axis];
%    Z_contour(1,1:2)=[Z_axis Z_axis];
%    contour_length_raw(1)=2;
%    psi_value_list_raw(1)=psi_scale(1);
%end
MIN_LENGTH=3

[min_val SMALLEST_CONTOUR]=min(contour_length_raw(2:end));
SMALLEST_CONTOUR=SMALLEST_CONTOUR+1

if SMALLEST_CONTOUR>1
    X_contour_raw(SMALLEST_CONTOUR+1,1:2)=[X_axis X_axis];
    Z_contour_raw(SMALLEST_CONTOUR+1,1:2)=[Z_axis Z_axis];
    contour_length_raw(SMALLEST_CONTOUR+1)=2;
    psi_value_list_raw(SMALLEST_CONTOUR+1)=psi_scale(1);
else
    psi_value_list_raw(2:Ncontours+1)=psi_value_list_raw(1:Ncontours);
    contour_length_raw(2:Ncontours+1)=contour_length_raw(1:Ncontours);
    X_contour_raw(2:Ncontours+1,:)=X_contour_raw(1:Ncontours,:);
    Z_contour_raw(2:Ncontours+1,:)=Z_contour_raw(1:Ncontours,:);
    X_contour_raw(1,1:2)=[X_axis X_axis];
    Z_contour_raw(1,1:2)=[Z_axis Z_axis];
    contour_length_raw(1)=2;
    psi_value_list_raw(1)=psi_scale(1);
end

index=1;

%if psi_value_list_raw(3)<psi_value_list_raw(end-1)
%    disp('reversing psi scale for contours.....')
%    psi_value_list_raw=flipud(psi_value_list_raw);
%	contour_length_raw=flipud(contour_length_raw);
%	X_contour_raw=flipud(X_contour_raw);
%	Z_contour_raw=flipud(Z_contour_raw);
%end



[min_val SMALLEST_CONTOUR]=min(contour_length_raw(1:end))


	

% psi_value_list=psi_value_list_raw(find(contour_length_raw>MIN_LENGTH));
% contour_length=contour_length_raw(find(contour_length_raw>MIN_LENGTH));






%for (n=2:4)
%index=index+1;
%for s=1:contour_length_raw(1)
%    X_contour(index,s)=X_contour_raw(index,s);
%    Z_contour(index,s)=Z_contour_raw(index,s);
%    contour_length(index)=contour_length_raw(index);
%    psi_value_list(index)=psi_value_list_raw(index);
%end
%end

if SIGN_CO_CURRENT_FIELD>0
    [psi_value_list_sorted sort_indexes]=sort(psi_value_list_raw,'ascend');
    X_contour_raw(:,:)=X_contour_raw(sort_indexes,:);
    Z_contour_raw(:,:)=Z_contour_raw(sort_indexes,:);
    contour_length_raw=contour_length_raw(sort_indexes);
    psi_value_list_raw=psi_value_list_sorted;
else
    [psi_value_list_sorted sort_indexes]=sort(psi_value_list_raw,'descend');
    X_contour_raw(:,:)=X_contour_raw(sort_indexes,:);
    Z_contour_raw(:,:)=Z_contour_raw(sort_indexes,:);
    contour_length_raw=contour_length_raw(sort_indexes);
    psi_value_list_raw=psi_value_list_sorted;
end


index=1;

%for s=1:contour_length_raw(1)
%    X_contour(index,s)=X_contour_raw(1,s);
%    Z_contour(index,s)=Z_contour_raw(1,s);
%    contour_length(index)=contour_length_raw(1);
%    psi_value_list(index)=psi_value_list_raw(1);
%end

for (n=1:length(contour_length_raw))
    if (contour_length_raw(n)>MIN_LENGTH) &&((psi_value_list_raw(n)~=psi_value_list_raw(n-1)))
        index=index+1;
        for s=1:contour_length_raw(n)
            X_contour(index,s)=X_contour_raw(n,s);
            Z_contour(index,s)=Z_contour_raw(n,s);
            contour_length(index)=contour_length_raw(n);
        end
        psi_value_list(index)=psi_value_list_raw(n);
    end
end
index


%%
[max_val LARGEST_CONTOUR]=max(contour_length(1:end))

if LARGEST_CONTOUR~=Nradial
%if LARGEST_CONTOUR<0.5*Nradial
%    NEW_LARGEST_CONTOUR=1
%else
	NEW_LARGEST_CONTOUR=Nradial
%end
    %Z_contour_raw_last=flipud(Z_contour_raw_last);
	%X_contour_raw_last=flipud(X_contour_raw_last);
    for s=1:min(length(X_contour_raw_last),length(Z_contour_raw_last))
        X_contour(NEW_LARGEST_CONTOUR,s)=X_contour_raw_last(s);
        Z_contour(NEW_LARGEST_CONTOUR,s)=Z_contour_raw_last(s);
    end

	contour_length(NEW_LARGEST_CONTOUR)=min(length(X_contour_raw_last),length(Z_contour_raw_last));
	psi_value_list(NEW_LARGEST_CONTOUR)=psi_scale(end);
	psi_value_list(NEW_LARGEST_CONTOUR-1)=0.5*(psi_scale(end)+psi_scale(end-1));
    % same as last contour, approximation for steep change of q
	for s=1:min(length(X_contour_raw_last),length(Z_contour_raw_last))
        X_contour(NEW_LARGEST_CONTOUR-1,s)=X_contour_raw_last(s);
        Z_contour(NEW_LARGEST_CONTOUR-1,s)=Z_contour_raw_last(s);
    end
	contour_length(NEW_LARGEST_CONTOUR-1)=min(length(X_contour_raw_last),length(Z_contour_raw_last));

	Ncontours=NEW_LARGEST_CONTOUR
end


if Ncontours~=Nradial
    disp('WARNING: the number of flux surfaces STILL does not match the number of contours !!!!');
end


%%

% now shifiting the data to separate the upper part 
% and the lower part of the flux surface contours

for (n=1:Ncontours)
    [residue index]=min(abs(Z_contour(n,:)-Z_axis));
    intersect_index(n)=index;
end

% the intersection is now the first index
X_contour_shift=X_contour;
Z_contour_shift=Z_contour;

for (n=2:Ncontours)
    index=contour_length(n)-intersect_index(n)+1;
    X_contour_shift(n,1:index)=X_contour(n,intersect_index(n):contour_length(n));
    Z_contour_shift(n,1:index)=Z_contour(n,intersect_index(n):contour_length(n));
    % extending the length to have an extra point for the intersection
    contour_length(n)=contour_length(n)+1;
    X_contour_shift(n,index+1:contour_length(n))=X_contour(n,1:intersect_index(n));
    Z_contour_shift(n,index+1:contour_length(n))=Z_contour(n,1:intersect_index(n));
end



% now looking for the middleintersection again
% but excluding the beginning and end of the contour
[residue index]=min(abs(Z_contour_shift(1,1:contour_length(n))-Z_axis));
intersect_index(1)=index;
for (n=2:2)
    [residue index]=min(abs(Z_contour_shift(n,3:contour_length(n)-1)-Z_axis));
    intersect_index(n)=index+2;
end
for (n=3:Ncontours)
    [residue index]=min(abs(Z_contour_shift(n,6:contour_length(n)-2)-Z_axis));
    intersect_index(n)=index+5;
end

X_contour_up=zeros(Ncontours,intersect_index(end));
Z_contour_up=zeros(Ncontours,intersect_index(end));

X_contour_down=zeros(Ncontours,intersect_index(end));
Z_contour_down=zeros(Ncontours,intersect_index(end));



for (n=2:Ncontours)
	if (Z_contour_shift(n,round(0.5*intersect_index(n)))<Z_axis)
		length_contour_down(n)=intersect_index(n);
		X_contour_down(n,1:intersect_index(n))=X_contour_shift(n,1:intersect_index(n));
		Z_contour_down(n,1:intersect_index(n))=Z_contour_shift(n,1:intersect_index(n));
	else
		length_contour_up(n)=intersect_index(n);
		X_contour_up(n,1:intersect_index(n))=X_contour_shift(n,1:intersect_index(n));
		Z_contour_up(n,1:intersect_index(n))=Z_contour_shift(n,1:intersect_index(n));
	end
    index=contour_length(n)-intersect_index(n);
    if (Z_contour_shift(n,round(0.5*intersect_index(n)+0.5*contour_length(n)))>Z_axis)
		length_contour_up(n)=index;
		X_contour_up(n,1:index)=X_contour_shift(n,intersect_index(n)+1:contour_length(n));
		Z_contour_up(n,1:index)=Z_contour_shift(n,intersect_index(n)+1:contour_length(n));
	else
		length_contour_down(n)=index;
		X_contour_down(n,1:index)=X_contour_shift(n,intersect_index(n)+1:contour_length(n));
		Z_contour_down(n,1:index)=Z_contour_shift(n,intersect_index(n)+1:contour_length(n));
	end
end
%for (n=2:Ncontours)
%    index=length_contour_up(n);
%    middle_up_index=round(0.5*index);
    % it may happen that the contour is not the upper part
%    if Z_contour_up(n,middle_up_index)<0
%        clear X_contour Z_contour
%        X_contour=X_contour_down;
%        Z_contour=Z_contour_down;
%        X_contour_down(n,1:index)=X_contour_up(n,1:index);
%        Z_contour_down(n,1:index)=Z_contour_up(n,1:index);
%        X_contour_up(n,1:length_contour_down(n))=X_contour(n,1:length_contour_down(n));
%        Z_contour_up(n,1:length_contour_down(n))=Z_contour(n,1:length_contour_down(n));
%        length_contour_up(n)=length_contour_down(n);
%        length_contour_down(n)=index;
%    end
%end

for (n=1:3)
    if length_contour_up(n)<=MIN_LENGTH
        X_contour_up(n,1:2)=[X_axis X_axis];
        Z_contour_up(n,1:2)=[Z_axis Z_axis];
        length_contour_up(n)=2;
    end
    if length_contour_down(n)<=MIN_LENGTH
        X_contour_down(n,1:2)=[X_axis X_axis];
        Z_contour_down(n,1:2)=[Z_axis Z_axis];
        length_contour_down(n)=2;
    end
end

max_contour_up=zeros(1,Ncontours);
min_contour_down=zeros(1,Ncontours);
max_contour_down=zeros(1,Ncontours);

for (n=2:Ncontours)
    [residue max_value]=max(X_contour_up(n,1:length_contour_up(n)));
    max_contour_up(n)=max_value;
    [residue min_value]=min(X_contour_up(n,1:length_contour_up(n)));
    min_contour_up(n)=min_value;
	if max_value<=min_value
		X_contour_up(n,1:length_contour_up(n))=fliplr(X_contour_up(n,1:length_contour_up(n)));
		Z_contour_up(n,1:length_contour_up(n))=fliplr(Z_contour_up(n,1:length_contour_up(n)));
	    [residue min_value]=min(X_contour_up(n,1:length_contour_up(n)));
		min_contour_up(n)=min_value;
		[residue max_value]=max(X_contour_up(n,1:length_contour_up(n)));
		max_contour_up(n)=max_value;
	end
    [residue min_value]=min(X_contour_down(n,1:length_contour_down(n)));
    min_contour_down(n)=min_value;
    [residue max_value]=max(X_contour_down(n,1:length_contour_down(n)));
    max_contour_down(n)=max_value;
	if max_value<=min_value
		X_contour_down(n,1:length_contour_down(n))=fliplr(X_contour_down(n,1:length_contour_down(n)));
		Z_contour_down(n,1:length_contour_down(n))=fliplr(Z_contour_down(n,1:length_contour_down(n)));
	    [residue min_value]=min(X_contour_down(n,1:length_contour_down(n)));
		min_contour_down(n)=min_value;
		[residue max_value]=max(X_contour_down(n,1:length_contour_down(n)));
		max_contour_down(n)=max_value;
	end
end




% small things to be adjusted on the extreme high field and low field
% side....


clear X_contour Z_contour
X_contour=X_contour_down;
Z_contour=Z_contour_down;

for (n=2:Ncontours)
    if min_contour_up(n)>1
        index=length_contour_down(n)+min_contour_up(n)-1;
        X_contour_down(n,length_contour_down(n)+1:index)=X_contour_up(n,1:min_contour_up(n)-1);
        Z_contour_down(n,length_contour_down(n)+1:index)=Z_contour_up(n,1:min_contour_up(n)-1);
        index=length_contour_up(n)-min_contour_up(n)+1;
        X_contour_up(n,1:index)=X_contour_up(n,min_contour_up(n):length_contour_up(n));
        Z_contour_up(n,1:index)=Z_contour_up(n,min_contour_up(n):length_contour_up(n));
        length_contour_up(n)=index;
    elseif max_contour_down(n)~=length_contour_down(n)
        index=length_contour_down(n)-max_contour_down(n)+1;
        X_contour_up(n,index+1:length_contour_up(n)+index)=X_contour_up(n,1:length_contour_up(n));
        Z_contour_up(n,index+1:length_contour_up(n)+index)=Z_contour_up(n,1:length_contour_up(n));
        X_contour_up(n,1:index)=X_contour_down(n,max_contour_down(n):length_contour_down(n));
        Z_contour_up(n,1:index)=Z_contour_down(n,max_contour_down(n):length_contour_down(n));
        length_contour_up(n)=length_contour_up(n)+index-1;
        length_contour_down(n)=length_contour_down(n)-index+1;
    end
end






% for (n=2:Ncontours)
%     [residue max_value]=max(X_contour_up(n,1:length_contour_up(n)));
%     max_contour_up(n)=max_value;
%     [residue min_value]=min(X_contour_down(n,1:length_contour_down(n)));
%     min_contour_down(n)=min_value;
%     [residue max_value]=max(X_contour_down(n,1:length_contour_down(n)));
%     max_contour_down(n)=max_value;
% end


% clear X_contour Z_contour
% X_contour=X_contour_up;
% Z_contour=Z_contour_up;
%
% for (n=2:Ncontours)
%     if max_contour_up(n)~=1
%        index=length_contour_up(n)+min_contour_down(n)-1;
%        X_contour_up(n,length_contour_up(n)+1:index)=X_contour_down(n,1:min_contour_down(n)-1);
%        Z_contour_up(n,length_contour_up(n)+1:index)=Z_contour_down(n,1:min_contour_down(n)-1);
%        index=length_contour_down(n)-min_contour_down(n)+1;
%        X_contour_down(n,1:index)=X_contour_down(n,min_contour_down(n):length_contour_down(n));
%        Z_contour_down(n,1:index)=Z_contour_down(n,min_contour_down(n):length_contour_down(n));
%        length_contour_down(n)=index;
%    end
% end

% clear X_contour Z_contour
% X_contour=X_contour_up;
% Z_contour=Z_contour_up;




% % removing identical Xvalues in the contour descriptions
% clear X_contour Z_contour
% X_contour=X_contour_up;
% Z_contour=Z_contour_up;
% 
% for (n=1:Ncontours)
%     u_contour=unique(X_contour_up(n,1:length_contour_up(n)));
%     val_contour=histc(X_contour_up(n,1:length_contour_up(n)),u_contour);
%     if max(val_contour)>1
%         if length(u_contour(val_contour>1))==1
%             duplicates_contour=find(X_contour_up(n,1:length_contour_up(n))==u_contour(val_contour>1));
%             del_index=duplicates_contour(1);
%             X_contour_up(n,1:length_contour_up(n)-1)=[X_contour(n,1:del_index-1) X_contour(n,del_index+1:length_contour_up(n))]; 
%             length_contour_up(n)=length_contour_up(n)-1;
%         end
%     end
% end

% checking for inverted X values on contour description
% for (n=2:Ncontours)
%     if X_contour_up(n,end)<X_contour_up(n,1)
%         X_contour_up(n,1:length_contour_up(n))=X_contour_up(n,length_contour_up(n):-1:1);
%         Z_contour_up(n,1:length_contour_up(n))=Z_contour_up(n,length_contour_up(n):-1:1);
%     end
%     if X_contour_down(n,end)<X_contour_down(n,1)
%         X_contour_down(n,1:length_contour_down(n))=X_contour_down(n,length_contour_down(n):-1:1);
%         Z_contour_down(n,1:length_contour_down(n))=Z_contour_down(n,length_contour_down(n):-1:1);
%     end
% end
 

% checking for non monotic contours X values at the end
% for (n=round(0.5*Ncontours):round(0.75*Ncontours))
%     cl=length_contour_up(n);
% 	s_inf=round(0.2*cl);
%     for s=s_inf:cl
%         if X_contour_up(n,s)<X_contour_up(n,s-1)
%             cl=min(cl,s-1);
%         end
%     end
%     length_contour_up(n)=cl;
% end
%  for (n=round(0.5*Ncontours):Ncontours-1)
%     cl=length_contour_down(n);
% 	s_inf=round(0.2*cl)+1;
%     for s=s_inf:cl
%         if X_contour_down(n,s)<X_contour_down(n,s-1)
%            cl=min(cl,s-1);
%         end
%     end
%     length_contour_down(n)=cl;
% end


for (n=1:Ncontours)
    for s=2:length_contour_down(n)-1
        if abs(X_contour_down(n,s)-X_contour_down(n,s-1))<1e-8
            X_contour_down(n,s)=0.5*(X_contour_down(n,s-1)+X_contour_down(n,s+1));
            Z_contour_down(n,s)=0.5*(Z_contour_down(n,s-1)+Z_contour_down(n,s+1));
        end
    end
    for s=2:length_contour_up(n)-1
        if abs(X_contour_up(n,s)-X_contour_up(n,s-1))<1e-8
            X_contour_up(n,s)=0.5*(X_contour_up(n,s-1)+X_contour_up(n,s+1));
            Z_contour_up(n,s)=0.5*(Z_contour_up(n,s-1)+Z_contour_up(n,s+1));
        end
    end
end
for (n=1:Ncontours)
    for s=2:length_contour_down(n)-1
        if abs(X_contour_down(n,s-1)-X_contour_down(n,s))<1e-8
            X_contour_down(n,s)=0.5*(X_contour_down(n,s+1)+X_contour_down(n,s-1));
            Z_contour_down(n,s)=0.5*(Z_contour_down(n,s+1)+Z_contour_down(n,s-1));
        end
    end
    for s=2:length_contour_up(n)-1
        if abs(X_contour_up(n,s-1)-X_contour_up(n,s))<1e-8
            X_contour_up(n,s)=0.5*(X_contour_up(n,s+1)+X_contour_up(n,s-1));
            Z_contour_up(n,s)=0.5*(Z_contour_up(n,s+1)+Z_contour_up(n,s-1));
        end
    end
end



% removing identical Xvalues in the contour descriptions
% clear X_contour Z_contour
% X_contour=X_contour_up;
% Z_contour=Z_contour_up;
% 
% for (n=1:Ncontours)
%     u_contour=unique(X_contour_up(n,1:length_contour_up(n)));
%     val_contour=histc(X_contour_up(n,1:length_contour_up(n)),u_contour);
%     if max(val_contour)>1
%         if length(u_contour(val_contour>1))==1
%             duplicates_contour=find(X_contour_up(n,1:length_contour_up(n))==u_contour(val_contour>1));
%             del_index=duplicates_contour(1);
%             X_contour_up(n,1:length_contour_up(n)-1)=[X_contour(n,1:del_index-1) X_contour(n,del_index+1:length_contour_up(n))]; 
%             length_contour_up(n)=length_contour_up(n)-1;
%         end
%     end
% end
 
hold on

for (n=1:Ncontours)
    plot(X_contour_down(n,1:length_contour_down(n)),Z_contour_down(n,1:length_contour_down(n)),'b');
    plot(X_contour_up(n,1:length_contour_up(n)),Z_contour_up(n,1:length_contour_up(n)),'r');
end

plot([X_axis X_axis],[Z_axis Z_axis],'+');
