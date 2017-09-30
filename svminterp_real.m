function h_FBMC_Aux = svminterp_real(index,hP_LS_FBMC_Aux,h,NrTime)
x=index.';
y=(hP_LS_FBMC_Aux).';
h_FBMC_Aux=nan(size(h));
model=svmtrain(y,x,'-s 3 -t 2 -c 2.2 -g 2.8 -p 0.01');
new_x=(1:NrTime).';
new_x(index)=[];
new_y=h(new_x);
[predict_real,mse_real,dec_real]=svmpredict(new_y,new_x,model);
h_FBMC_Aux(new_x)=predict_real;
h_FBMC_Aux(index)=hP_LS_FBMC_Aux;
end