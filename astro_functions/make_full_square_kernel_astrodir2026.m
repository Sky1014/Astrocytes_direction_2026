function Full_Square_Kernel = make_full_square_kernel_astrodir2026(stage, tr)
    square_wave_stage = [];
    square_wave_stage = ones(tr(3,end),1);
    for i = 1:length(stage(1,:))
        square_wave_stage(stage(2,i):stage(3,i),1) = 1.05;
    end
    Full_Square_Kernel = square_wave_stage;
end
