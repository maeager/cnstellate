

for ii = 0:5
 vs = load (["vsSPIKES." num2str(ii) ".dat"]);
 ray = load(["rayltest." num2str(ii) ".dat"]);
 significant = ray > 5.06;
 z0 = significant .* vs;
 subplot(3,2,ii+1)
 imagesc([0:99],0:5:70,z0', [0 1]);axis("xy"); 

end
xlim([50 100]); 

 xlim([50 100]); 
 xlabel (" Channel No.");
 ylabel (" Level (dB SPL)" ); 
