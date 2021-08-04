function [mask_overlap,mask_union,mask_overlap_union, SO_detected] = SO_overlap_detect(s1,s2,fs,SO_length)
% mask_overlap: 
% s1: channel x time
[dim,~]=size(s1);
mask_overlap= zeros(size(s1));
mask_union = zeros(size(s1));
SO_detected = cell(dim,2);

for i=1:dim
    [SO_mask1,SO_de1] = SO_detect(s1(i,:),[0.5,5],0.5,fs);
    [SO_mask2,SO_de2] = SO_detect(s2(i,:),[0.5,5],0.5,fs);
    mask_overlap(i,:) =SO_mask1.*SO_mask2; 
    mask_union(i,:) =double(SO_mask1|SO_mask2);
    
    mask_overlap_union = zeros(size(mask_union));
    st = 0;
    ed = 0;
    for kk=1:length(mask_union)-1
        if(mask_union(1,kk)==0)&&(mask_union(1,kk+1)==1)
            st = kk;
        end
        if(mask_union(1,kk)==1)&&(mask_union(1,kk+1)==0)
            ed = kk;
        end
        if st<ed
            if sum(mask_overlap(1,st:ed))>0
                mask_overlap_union(1,st:ed)=1;
            end
        end
    end
    
    
    SO_1 =[];
    SO_2 = [];
    st = 0;
    ed=0;
    flag = 0;
    for kk=1:length(mask_overlap_union)-1
        if(mask_overlap_union(1,kk)==0)&&(mask_union(1,kk+1)==1)
            st = kk;
            flag=1;
        end
        if(mask_union(1,kk)==1)&&(mask_union(1,kk+1)==0)
            ed = kk;
            flag=2;
        end
        if flag==2
            t_middle = round((st+ed)/2);
            if t_middle-SO_length>=1&&t_middle+SO_length<=length(s1)
                SO_1_tmp = s1(i,t_middle-SO_length:t_middle+SO_length);
                SO_2_tmp = s2(i,t_middle-SO_length:t_middle+SO_length);
                SO_1 = [SO_1;SO_1_tmp];
                SO_2 = [SO_2;SO_2_tmp];
            end
        end
        flag=0;
    end
    SO_detected{i,1}=SO_1;
    SO_detected{i,2}=SO_2;
end




