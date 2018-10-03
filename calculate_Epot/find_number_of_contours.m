function nb_contours=find_number_of_contours(Contours_object)

index=1;
nb_contours=0;
index_max=size(Contours_object,2);
while index<index_max
    nb_contours=nb_contours+1;
    index=index+Contours_object(2,index)+1;
end

return