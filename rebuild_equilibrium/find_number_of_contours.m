function nb_contours=find_number_of_contours(Contours_object)

index=1;
nb_contours=0;
index_max=size(Contours_object,2);
while index<index_max
    size_contour=Contours_object(2,index);
    if size_contour>6
        nb_contours=nb_contours+1;
    end
    index=index+size_contour+1;
end

return