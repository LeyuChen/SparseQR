
load covariate_specification.mat;

s1="(";s2=",";s3=")";

if p==21
vname=v21;

elseif p==49
vname = v49;

else

vname = string(zeros(p,1));
B_x = string(zeros(4*B_terms,1)); % mage, fage, nprenatal, monthslb
for k=1:4
    switch k
        case 1
            spline_var = "mage";
        case 2
            spline_var = "fage";
        case 3
            spline_var = "npre";
        case 4
            spline_var = "mon";
    end
    for j=1:B_terms
     B_x(B_terms*(k-1)+j,1)=strcat("B",num2str(j),s1,spline_var,s3);  
    end
end

vname(1:21,1)=table2array(v49(1:21,1)); vname(22:21+length(B_x),1)=B_x;
ind = 21+length(B_x);
for j=1:20
 for k=1:length(B_x)
 vname(ind+(j-1)*length(B_x)+k)=strcat(vname(j+1),"*",B_x(k));
 end  
end    

end

sel_L=(abs(bhat_L)>tol);sel_U=(abs(bhat_U)>tol);

for i=1:p
   bsel_L(i) = mean(abs(bhat_L(i,:))>tol); 
   bsel_U(i) = mean(abs(bhat_U(i,:))>tol); 
end

for i=1:p
coef_L(i,1) = mean(bhat_L(i,sel_L(i,:))); 
coef_U(i,1) = mean(bhat_U(i,sel_U(i,:)));
end

tlist = 1:10;
[sorted_L,ind_L]=sort(bsel_L,'descend');
[sorted_U,ind_U]=sort(bsel_U,'descend');



if p<=49
v_L=table2array(vname(ind_L(tlist),1));v_U=table2array(vname(ind_U(tlist),1));
else
v_L=vname(ind_L(tlist),1);v_U=vname(ind_U(tlist),1);
end

prop_L = string(round(sorted_L(tlist),2)); prop_U = string(round(sorted_U(tlist),2)); 
sorted_coef_L = string(round(coef_L(ind_L(tlist)),2)); sorted_coef_U = string(round(coef_U(ind_U(tlist)),2));
r_L=strcat(s1,prop_L,s2,sorted_coef_L,s3); r_U=strcat(s1,prop_U,s2,sorted_coef_U,s3);
result_L = string(zeros(2*length(tlist),1));result_U=string(zeros(2*length(tlist),1));
result_L(1:2:19)=v_L;result_L(2:2:20)=r_L;
result_U(1:2:19)=v_U;result_U(2:2:20)=r_U;

disp(['Variable selection Results for ' QR_method ' with p = ' num2str(p)]);
disp('Top 10 most often selected variables at quantile level 0.05');
disp(result_L);
disp('Top 10 most often selected variables at quantile level 0.95');
disp(result_U);
