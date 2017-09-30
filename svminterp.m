function h_FBMC_Aux = svminterp(index,hP_LS_FBMC_Aux,h,NrTime)
x=index.';
y_magni=abs(hP_LS_FBMC_Aux).';
y_phase=angle(hP_LS_FBMC_Aux).';
h_FBMC_Aux=nan(size(h));
model_magni=svmtrain(y_magni,x,'-s 3 -t 2 -c 2.2 -g 2.8 -p 0.01');
model_phase=svmtrain(y_phase,x,'-s 3 -t 2 -c 2.2 -g 2.8 -p 0.01');
new_x=(1:NrTime).';
new_x(index)=[];
new_magni=abs(h(new_x));
new_phase=angle(h(new_x));
[predict_magni,mse_real,dec_real]=svmpredict(new_magni,new_x,model_magni);
[predict_phase,mse_img,dec_img]=svmpredict(new_phase,new_x,model_phase);
h_FBMC_Aux(new_x)=predict_magni.*exp(1j*predict_phase);
h_FBMC_Aux(index)=hP_LS_FBMC_Aux;
end