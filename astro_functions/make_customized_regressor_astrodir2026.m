function Customized_Regressor = make_customized_regressor_astrodir2026(k, s, f, c, r)
    regressor = k;
    Regressor = zeros(1,length(f));
    Regressor = regressor(f);
    Rho = corr(Regressor, c.');
    Customized_Regressor.regressor = Regressor;
    if r >= 0
        Customized_Regressor.cell_index = find(Rho > r);
        Customized_Regressor.rho = Rho(find(Rho > r));
    else
        Customized_Regressor.cell_index = find(Rho <= r);
        Customized_Regressor.rho = Rho(find(Rho <= r));
    end

    cr = c(Customized_Regressor.cell_index,:);
    NA_Rows = find(any(isnan(cr),2) == 1);
    Customized_Regressor.cell_index(NA_Rows) = [];
    Customized_Regressor.rho(NA_Rows) = [];

end
