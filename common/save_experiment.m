function method_list = save_experiment(method_list, i, index_fail, num)
    % save result
    for k= 1:length(method_list)
       %results without deleting solutions
        tmethod_list = method_list(k);
        method_list(k).c(index_fail{k}) = [];
        method_list(k).e(index_fail{k}) = [];
        method_list(k).r(index_fail{k}) = [];
        method_list(k).t(index_fail{k}) = [];

        method_list(k).pfail(i) = 100 * numel(index_fail{k})/num;

        method_list(k).mean_c(i)= mean(method_list(k).c);
        method_list(k).mean_e(i)= mean(method_list(k).e);
        method_list(k).med_c(i)= median(method_list(k).c);
        method_list(k).med_e(i)= median(method_list(k).e);
        method_list(k).std_c(i)= std(method_list(k).c);
        method_list(k).std_e(i)= std(method_list(k).e);

        method_list(k).mean_r(i)= mean(method_list(k).r);
        method_list(k).mean_t(i)= mean(method_list(k).t);
        method_list(k).med_r(i)= median(method_list(k).r);
        method_list(k).med_t(i)= median(method_list(k).t);
        method_list(k).std_r(i)= std(method_list(k).r);
        method_list(k).std_t(i)= std(method_list(k).t);

        %results deleting solutions where not all the methods finds one
        tmethod_list.c(unique([index_fail{:}])) = [];
        tmethod_list.e(unique([index_fail{:}])) = [];
        tmethod_list.r(unique([index_fail{:}])) = [];
        tmethod_list.t(unique([index_fail{:}])) = [];

        method_list(k).deleted_mean_c(i)= mean(tmethod_list.c);
        method_list(k).deleted_mean_e(i)= mean(tmethod_list.e);
        method_list(k).deleted_med_c(i)= median(tmethod_list.c);
        method_list(k).deleted_med_e(i)= median(tmethod_list.e);
        method_list(k).deleted_std_c(i)= std(tmethod_list.c);
        method_list(k).deleted_std_e(i)= std(tmethod_list.e);

        method_list(k).deleted_mean_r(i)= mean(tmethod_list.r);
        method_list(k).deleted_mean_t(i)= mean(tmethod_list.t);
        method_list(k).deleted_med_r(i)= median(tmethod_list.r);
        method_list(k).deleted_med_t(i)= median(tmethod_list.t);
        method_list(k).deleted_std_r(i)= std(tmethod_list.r);
        method_list(k).deleted_std_t(i)= std(tmethod_list.t);
    end
end

