function result = b_fil(pos)



     result = pos;
     A = -0.9629940509502;
     B = [0.018502974524892, 0.018502974524892]; 
     num_col = size(pos, 2); 
     num_row = size(pos, 1);
     for i = 1:num_row
	    for j= 1:num_col
	        if (j == 1) 
		        result(i,j) = B(1) * pos(i,j) + B(2) * 0 - A*0;
            else
		        result(i,j) = B(1) * pos(i,j) + B(2) * pos(i, j-1) - A* result(i,j-1);
            end
        end
     end
     
end