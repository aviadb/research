addpath('/home/aviadb/Desktop/thesis/sap-voicebox/voicebox/');

clc;
close all;





% f = 1:100:22e3;
% wav_dir = '/home/aviadb/Desktop/study/525.804_Thesis_II/research/Scripts/enhanced_metrics/';
% close_mic = audioread([wav_dir, 'F05_440C020E_BUS.CH0.wav']);
% noisy = audioread([wav_dir, 'F05_440C020E_BUS.CH1.wav']);
% enhanced = audioread([wav_dir, 'F05_440C020E_BUS.wav']);
% enhanced2 = audioread([wav_dir, 'F05_440C020E_BUS_enh2.wav']);

% enhanced_scaled = (0.5-(-0/5)) * (enhanced - min(enhanced))/(max(enhanced) - min(enhanced)) + (-0.5);
% snrseg(noisy, close_mic, 16e3, 'wp');

% figure();
% plot(enhanced);
% hold on;
% plot(close_mic + 0.2);
close all;
% figure();
% v_sigalign(enhanced, close_mic, 100, 'uqp', 16000);

% figure();
% v_snrseg(enhanced, close_mic, 16e3, 'Vqp');
% 
% figure();
% v_snrseg(enhanced./2, close_mic, 16e3, 'wqp');

% figure();
% [seg_enh, glob_enh] = v_snrseg(enhanced./2, close_mic, 16e3, 'vp');
% 
% 
% % figure();
% % v_snrseg(enhanced_scaled, close_mic, 16e3, 'wqp');
% 
% 
% figure();
% [seg_noisy, glob_noisy] = v_snrseg(noisy, close_mic, 16e3, 'Vp');
% % 
% % % figure();
% % snrseg(enhanced, close_mic(84:end), 16e3, 'Vp');
% % 
% % figure();
% % snrseg(enhanced, close_mic, 16e3, 'Vp');
% % 
% % figure();
% % v_sigalign(enhanced./2, close_mic, 100, 'up', 16000);
% 
% % figure();
% % v_snrseg(enhanced./2, noisy, 16e3, 'Vqp');
% 
% seg_snr_improve = abs(seg_noisy-seg_enh)
% snr_improve = abs(glob_noisy-glob_enh)


% wav_dir = '/home/aviadb/Desktop/thesis/beamform_NN_IRM/reco_enh/F04_050C0109_BUS/';
% close_mic = audioread([wav_dir, 'F04_050C0109_BUS.CH2.fclean.wav']);
% noisy = audioread([wav_dir, 'F04_050C0109_BUS.CH2.wav']);
% ideal_irm = audioread([wav_dir, 'F04_050C0109_BUS_IRM_ch5_enh.wav']);



% nn_simu_enh_dir = '/home/aviadb/Desktop/thesis/beamform_NN_IRM/enhanced/NN_IRM/dt05_bus_simu/';
% nn_real_enh_dir = '/home/aviadb/Desktop/thesis/beamform_NN_IRM/enhanced/NN_IRM/dt05_bus_real/';

% F06_447C0217
% F04_050C0109



% close all;


mix_dir = '/home/aviadb/Desktop/datasets/CHiME3/data/audio/16kHz/isolated/';
s_clean_dir = '/home/aviadb/Desktop/datasets/CHiME3/data/audio/16kHz/isolated_fn/';
real_dir = '/home/aviadb/Desktop/datasets/CHiME3/data/audio/16kHz/isolated/';
enh_dir = '/home/aviadb/Desktop/thesis/beamform_NN_IRM/enhanced/NN_IRM/';

snr_df = table('Size',[0 34], ...
    'VariableNames',{'seg_enh', 'seg_noisy', 'glob_noisy', 'glob_enh', ...
                    'SNR', 'Seg_SNR','SNR_sign', 'Seg_SNR_sign',...
                    'mode', 'channel', 'enh', 'set', 'scenario', 'utt',...
                    'Ideal_IRM_Seg_SNR', 'Ideal_IRM_SNR',...
                    'Ideal_cIRM_Seg_SNR', 'Ideal_cIRM_SNR',...
                    'Ideal_PSM_Seg_SNR', 'Ideal_PSM_SNR',...
                    'Ideal_ORM_Seg_SNR', 'Ideal_ORM_SNR',...
                    'cIRM_SNR',...
                    'cIRM_Seg_SNR',...
                    'PSM_SNR',...
                    'PSM_Seg_SNR',...
                    'ORM_SNR',...
                    'ORM_Seg_SNR',...
                    'cIRM_SNR_Ratio', 'cIRM_Seg_SNR_Ratio',...
                    'PSM_SNR_Ratio', 'PSM_Seg_SNR_Ratio',...
                    'ORM_SNR_Ratio', 'ORM_Seg_SNR_Ratio',...
                    }, ...
    'VariableTypes',{'double','double','double','double',...
                    'double','double','double','double',...
                    'string', 'uint8', 'string', 'string', 'string', 'string',...
                    'double','double',...
                    'double','double',...
                    'double','double',...
                    'double','double',...
                    'double',...
                    'double',...
                    'double',...
                    'double',...
                    'double',...
                    'double',...
                    'double','double',...
                    'double','double',...
                    'double','double',...
                    });
plot_flag = false;

scenarios = {'BUS'};
sets = {'dt', 'et'};

dt_utts = {'F04_050C0109', 'M04_423C0216', 'F01_423C0213', 'M03_423C0211'};
et_utts = {'F06_447C0217', 'M06_447C0216', 'M05_447C0215', 'F05_447C0211'};
% utts = {'F04_050C0109', 'F06_447C0217'};
modes = {'real', 'simu'};

num_channels = 6;
frame_length = 400;
fft_length   = 512;
frame_shift  = 100;
hanning_wnd = v_windows(3,frame_length,'l')';

for scn_idx=1:length(scenarios)
    scen = scenarios{scn_idx};
    
    for set_idx=1:length(sets)
        set = sets{set_idx};
        
        for mode_idx=1:length(modes)
            mode = modes{mode_idx};
            dir_pre = [set,'05_',lower(scen),'_',mode,'/'];

            utts =  dir([enh_dir dir_pre '*.wav']);
            utts =  vertcat(utts.name);
%             for utt_idx=1:length(utts)
            for utt_idx=1:10 % debug
                if mod(utt_idx, 50) == 0
                    sprintf('set: %d(%d), mode: %d(%d), utt: %d(%d)', ...
                        set_idx, length(sets), mode_idx, ...
                        length(modes), utt_idx, length(utts))
                end
%                 utt = [utts{utt_idx} '_' scen];
                utt = utts(utt_idx,1:end-4);

                enh_irm = audioread([enh_dir, dir_pre, utt, '.wav']);
    
    %             s_enh_irm = audioread([enh_dir, utt, '.wav']);
        
    %             r_enh_irm = audioread([enh_dir, dir_pre 'real/' F04_050C0109_BUS.wav']);
    %             r_clean = audioread([real_dir, 'F04_050C0109_BUS.CH0.wav']);
                if mode=='real'
                    clean = audioread([real_dir, dir_pre, utt, '.CH0.wav']);
                end
%                 for ch=1:6
                ch = 5; % Only channel 5
                    if mode=='simu'
                        clean = audioread([s_clean_dir,dir_pre,utt,'.CH',num2str(ch),'.fclean.wav']);
                    end
    
                    noisy = audioread([mix_dir,dir_pre,utt,'.CH',num2str(ch),'.wav']);
                    % Real Wavs
    %                 r_nsy_filename = [
    %                     dir_pre,'_real/',utt,'.CH',num2str(i),'.wav'
    %                 ];
    %                 r_nsy_filename = ['F04_050C0109_BUS.CH' num2str(i) '.wav'];
    %                 r_nsy = audioread([real_dir, r_nsy_filename]);
    %             
                    % Simu Wavs
    %                 s_clean_filename = ['F04_050C0109_BUS.CH' num2str(i) '.fclean.wav'];
    %                 s_clean = audioread([s_clean_dir, s_clean_filename]);
    %                 s_nsy_filename = ['F04_050C0109_BUS.CH' num2str(i) '.wav'];
    %                 s_nsy = audioread([simu_dir, s_nsy_filename]);
    
                    if (plot_flag)
                        figure();
                        [seg_enh, glob_enh] = v_snrseg(enh_irm, clean, 16e3, 'vp');
                        figure();
                        [seg_noisy, glob_noisy] = v_snrseg(noisy, clean, 16e3, 'Vp');
                    else
                        [seg_enh, glob_enh] = v_snrseg(enh_irm, clean, 16e3, 'v');
                        [seg_noisy, glob_noisy] = v_snrseg(noisy, clean, 16e3, 'V');
                    end
    
    %                 if (plot_flag)
    %                     figure();
    %                     [r_seg_enh, r_glob_enh] = v_snrseg(r_enh_irm, r_clean, 16e3, 'vp');
    %                     figure();
    %                     [r_seg_noisy, r_glob_noisy] = v_snrseg(r_nsy, r_clean, 16e3, 'Vp');
    %             
    %                     figure();
    %                     [s_seg_enh, s_glob_enh] = v_snrseg(s_enh_irm, s_clean, 16e3, 'vp');
    %                     figure();
    %                     [s_seg_noisy, s_glob_noisy] = v_snrseg(s_nsy, s_clean, 16e3, 'Vp');
    %                 else
    %                     [r_seg_enh, r_glob_enh] = v_snrseg(r_enh_irm, r_clean, 16e3, 'v');
    %                     [r_seg_noisy, r_glob_noisy] = v_snrseg(r_nsy, r_clean, 16e3, 'V');
    %             
    %                     [s_seg_enh, s_glob_enh] = v_snrseg(s_enh_irm, s_clean, 16e3, 'v');
    %                     [s_seg_noisy, s_glob_noisy] = v_snrseg(s_nsy, s_clean, 16e3, 'V');
    %                 end
                %     if ~sign(glob_noisy)
                %         glob_noisy = -glob_noisy;
                %     end
                    seg_snr_imp = abs(seg_noisy-seg_enh);
                    if seg_noisy > seg_enh
                        seg_snr_imp_sign = -seg_snr_imp;
                    else
                        seg_snr_imp_sign = seg_snr_imp;
                    end
                    snr_imp = abs(glob_noisy-glob_enh);
                    if glob_noisy > glob_enh
                        snr_imp_sign = -snr_imp;
                    else
                        snr_imp_sign = snr_imp;
                    end


                    clean = audioread([s_clean_dir,dir_pre,utt,'.CH',num2str(ch),'.fclean.wav']);
                    noise = audioread([s_clean_dir dir_pre utt '.CH' int2str(ch) '.noise.wav']);
    

                    % Extract ideals
                    noisy_stft = stft_voicebox(noisy);
                    clean_stft = stft_voicebox(clean);
                    if mode=='real'%real case
                        noise = estimate_noise_v2(clean, noisy);
                        noise_stft = stft_voicebox(noise);
                    else
                        noise_stft = stft_voicebox(noise);
                    end

                    
                    alpha = 0.5;
                    ideal_IRM_clean = ((abs(clean_stft).^2)./(abs(clean_stft).^2+abs(noise_stft).^2)).^alpha;
                    ideal_IRM_noise = ((abs(noise_stft).^2)./(abs(clean_stft).^2+abs(noise_stft).^2)).^alpha;
                
                    Yr = real(noisy_stft);
                    Yi = imag(noisy_stft);
                    Sr = real(clean_stft);
                    Si = imag(clean_stft);
                    Nr = real(noise_stft);
                    Ni = imag(noise_stft);
                    ideal_cIRM_clean_r = ((Yr.*Sr) + (Yi.*Si))./(Yr.^2 + Yi.^2);
                    ideal_cIRM_clean_i = ((Yr.*Si) - (Yi.*Sr))./(Yr.^2 + Yi.^2);
                    ideal_cIRM_clean = ideal_cIRM_clean_r + 1j*ideal_cIRM_clean_i;
                
                    ideal_cIRM_noise_r = ((Yr.*Nr) + (Yi.*Ni))./(Yr.^2 + Yi.^2);
                    ideal_cIRM_noise_i = ((Yr.*Ni) + (Yi.*Nr))./(Yr.^2 + Yi.^2);
                    ideal_cIRM_noise = ideal_cIRM_noise_r + 1j*ideal_cIRM_noise_i;
                    
                    ideal_PSM_clean = (abs(clean_stft) ./ abs(noisy_stft)) .* cos(angle(clean_stft) - angle(noisy_stft));
                    ideal_PSM_noise = (abs(noise_stft) ./ abs(noisy_stft)) .* cos(angle(noise_stft) - angle(noisy_stft));
                        
                    ideal_ORM_clean = ((abs(clean_stft).^2) + real(clean_stft .* conj(noise_stft)));
                    ideal_ORM_clean = ideal_ORM_clean ./ ((abs(clean_stft).^2) + (abs(noise_stft).^2) + 2*real(clean_stft .* conj(noise_stft)));
                    ideal_ORM_noise = ((abs(noise_stft).^2) + real(noise_stft .* conj(clean_stft)));
                    ideal_ORM_noise = ideal_ORM_noise ./ ((abs(clean_stft).^2) + (abs(noise_stft).^2) + 2*real(noise_stft .* conj(clean_stft)));
                    
                    % Multiply and Reconstruct
                    ideal_irm_S = ideal_IRM_clean'.*noisy_stft'; 
                    frames_enhan = irfft(ideal_irm_S, fft_length, 2);
                    ideal_irm_enh = overlapadd(frames_enhan(:, 1: frame_length), hanning_wnd, frame_shift);

                    [id_irm_seg_enh, id_irm_glob_enh] = v_snrseg(ideal_irm_enh, clean, 16e3, 'v');
                    
                    
                    ideal_cirm_S = ideal_cIRM_clean'.*noisy_stft'; 
                    frames_enhan = irfft(ideal_cirm_S, fft_length, 2);
                    ideal_cirm_enh = overlapadd(frames_enhan(:, 1: frame_length), hanning_wnd, frame_shift);

                    [id_cirm_seg_enh, id_cirm_glob_enh] = v_snrseg(ideal_cirm_enh, clean, 16e3, 'v');
                       

                    ideal_psm_S = ideal_PSM_clean'.*noisy_stft'; 
                    frames_enhan = irfft(ideal_psm_S, fft_length, 2);
                    ideal_psm_enh = overlapadd(frames_enhan(:, 1: frame_length), hanning_wnd, frame_shift);

                    [id_psm_seg_enh, id_psm_glob_enh] = v_snrseg(ideal_psm_enh, clean, 16e3, 'v');
                       

                    ideal_orm_S = ideal_ORM_clean'.*noisy_stft'; 
                    frames_enhan = irfft(ideal_orm_S, fft_length, 2);
                    ideal_orm_enh = overlapadd(frames_enhan(:, 1: frame_length), hanning_wnd, frame_shift);

                    [id_orm_seg_enh, id_orm_glob_enh] = v_snrseg(ideal_orm_enh, clean, 16e3, 'v');
                       
                    [seg_noisy, glob_noisy] = v_snrseg(noisy, clean, 16e3, 'V');


                    id_irm_seg_snr_imp = abs(seg_noisy-id_irm_seg_enh);
%                     if seg_noisy > seg_enh
%                         id_irm_seg_snr_imp_sign = -seg_snr_imp;
%                     else
%                         id_irm_seg_snr_imp_sign = seg_snr_imp;
%                     end
                    id_irm_snr_imp = abs(glob_noisy-id_irm_glob_enh);
%                     if glob_noisy > glob_enh
%                         id_irm_snr_imp_sign = -snr_imp;
%                     else
%                         id_irm_snr_imp_sign = snr_imp;
%                     end

                    id_cirm_seg_snr_imp = abs(seg_noisy-id_cirm_seg_enh);
                    id_cirm_snr_imp = abs(glob_noisy-id_cirm_glob_enh);

                    id_psm_seg_snr_imp = abs(seg_noisy-id_psm_seg_enh);
                    id_psm_snr_imp = abs(glob_noisy-id_psm_glob_enh);

                    id_orm_seg_snr_imp = abs(seg_noisy-id_orm_seg_enh);
                    id_orm_snr_imp = abs(glob_noisy-id_orm_glob_enh);                    
    %                 s_seg_snr_imp = abs(s_seg_noisy-s_seg_enh);
    %                 s_snr_imp = abs(s_glob_noisy-s_glob_enh);

%                     cellData = {
%                         seg_enh, seg_noisy, glob_noisy, glob_enh, ...
%                         snr_imp, seg_snr_imp, snr_imp_sign, seg_snr_imp_sign,...
%                          mode, ch, 'IRM', set, 'BUS', utt,...
%                          id_irm_seg_snr_imp, id_irm_snr_imp,...
%                          id_cirm_seg_snr_imp, id_cirm_snr_imp,...
%                          id_psm_seg_snr_imp, id_psm_snr_imp,...
%                          id_orm_seg_snr_imp, id_orm_snr_imp,...
%                          snr_imp * ((id_cirm_snr_imp/id_irm_snr_imp*100)-100)/100,...
%                          seg_snr_imp * ((id_cirm_seg_snr_imp/id_irm_seg_snr_imp*100)-100)/100,...                         
%                          snr_imp * ((id_psm_snr_imp/id_irm_snr_imp*100)-100)/100,...
%                          seg_snr_imp * ((id_psm_seg_snr_imp/id_irm_seg_snr_imp*100)-100)/100,...
%                          snr_imp * ((id_orm_snr_imp/id_irm_snr_imp*100)-100)/100,...
%                          seg_snr_imp * ((id_orm_seg_snr_imp/id_irm_seg_snr_imp*100)-100)/100,...
%                          ((id_cirm_snr_imp/id_irm_snr_imp*100)-100)/100, ...
%                          ((id_cirm_seg_snr_imp/id_irm_seg_snr_imp*100)-100)/100, ...
%                          ((id_psm_snr_imp/id_irm_snr_imp*100)-100)/100,...
%                          ((id_psm_seg_snr_imp/id_irm_seg_snr_imp*100)-100)/100,...
%                          ((id_orm_snr_imp/id_irm_snr_imp*100)-100)/100,...
%                          ((id_orm_seg_snr_imp/id_irm_seg_snr_imp*100)-100)/100,...
%                          };
%                     snr_df = [snr_df;cellData];
%                 end % End for channel
            end
        end
    end
end
mean(snr_df.SNR(snr_df.mode=='real' & snr_df.set=='dt'))
mean(snr_df.SNR_sign(snr_df.mode=='real' & snr_df.set=='dt'))
mean(snr_df.Seg_SNR(snr_df.mode=='real' & snr_df.set=='dt'))
mean(snr_df.Seg_SNR_sign(snr_df.mode=='real' & snr_df.set=='dt'))

r_dt_seg_snr_t = [snr_df.Seg_SNR(snr_df.mode=='real' & snr_df.set=='dt'), snr_df.Seg_SNR_sign(snr_df.mode=='real' & snr_df.set=='dt')];
r_dt_snr_t = [snr_df.SNR(snr_df.mode=='real' & snr_df.set=='dt'), snr_df.SNR_sign(snr_df.mode=='real' & snr_df.set=='dt')];

r_et_seg_snr_t = [snr_df.Seg_SNR(snr_df.mode=='real' & snr_df.set=='et'), snr_df.Seg_SNR_sign(snr_df.mode=='real' & snr_df.set=='et')];
r_et_snr_t = [snr_df.SNR(snr_df.mode=='real' & snr_df.set=='et'), snr_df.SNR_sign(snr_df.mode=='real' & snr_df.set=='et')];

s_dt_seg_snr_t = [snr_df.Seg_SNR(snr_df.mode=='simu' & snr_df.set=='dt'), snr_df.Seg_SNR_sign(snr_df.mode=='simu' & snr_df.set=='dt')];
s_dt_snr_t = [snr_df.SNR(snr_df.mode=='simu' & snr_df.set=='dt'), snr_df.SNR_sign(snr_df.mode=='simu' & snr_df.set=='dt')];

s_et_seg_snr_t = [snr_df.Seg_SNR(snr_df.mode=='simu' & snr_df.set=='et'), snr_df.Seg_SNR_sign(snr_df.mode=='simu' & snr_df.set=='et')];
s_et_snr_t = [snr_df.SNR(snr_df.mode=='simu' & snr_df.set=='et'), snr_df.SNR_sign(snr_df.mode=='simu' & snr_df.set=='et')];


mean(r_dt_seg_snr_t)
mean(r_dt_seg_snr_t(:))

mean(r_dt_snr_t)
mean(r_dt_snr_t(:))

mean(r_et_seg_snr_t)
mean(r_et_seg_snr_t(:))

mean(r_et_snr_t)
mean(r_et_snr_t(:))



mean(s_dt_seg_snr_t)
mean(s_dt_seg_snr_t(:))

mean(s_dt_snr_t)
mean(s_dt_snr_t(:))

mean(s_et_seg_snr_t)
mean(s_et_seg_snr_t(:))

mean(s_et_snr_t)
mean(s_et_snr_t(:))

df_save_dir = '/home/aviadb/Desktop/thesis/speechbrain/recipes/CommonVoice/ASR/transformer/';
parquetwrite([df_save_dir 'snr_segsnr_df.parquet'], snr_df);


%     s_noisy_filename = ['F04_050C0109_BUS.CH' num2str(i) '.wav'];
%     s_noisy = audioread([simu_dir, s_noisy_filename]);
%     clean_filename = ['F04_050C0109_BUS.CH' num2str(i) '.fclean.wav'];
%     close_mic = audioread([wav_dir, clean_filename]);
%     figure();
%     [simu_seg_enh, simu_glob_enh] = v_snrseg(enh_irm, close_mic, 16e3, 'v');
%     figure();
%     [simu_seg_noisy, simu_glob_noisy] = v_snrseg(noisy, close_mic, 16e3, 'V');
%     if ~sign(glob_noisy)
%         glob_noisy = -glob_noisy;
%     end




% for i=1:6
%     noisy_filename = ['F04_050C0109_BUS.CH' num2str(i) '.wav'];
%     noisy = audioread([simu_dir, noisy_fiflename]);
%     figure();
%     [seg_enh, glob_enh] = v_snrseg(enh_irm, close_mic, 16e3, 'vp');
%     figure();
%     [seg_noisy, glob_noisy] = v_snrseg(noisy, close_mic, 16e3, 'Vp');
%     if ~sign(glob_noisy)
%         glob_noisy = -glob_noisy;
%     end
%     seg_snr_improve = abs(seg_noisy-seg_enh)
%     snr_improve = abs(glob_noisy-glob_enh)
%     cellData = {snr_improve, seg_snr_improve,'real',i,'IRM'};
%     snr_df = [snr_df;cellData];
% end

% figure();
% [seg_enh, glob_enh] = v_snrseg(ideal_irm, close_mic, 16e3, 'vp');
% figure();
% [seg_noisy, glob_noisy] = v_snrseg(noisy, close_mic, 16e3, 'Vp');
% seg_snr_improve = abs(seg_noisy-seg_enh)
% snr_improve = abs(glob_noisy-glob_enh)
% 
% figure();
% [seg_enh, glob_enh] = v_snrseg(enh_irm, close_mic, 16e3, 'qwp', 0.00625);
% figure();
% [seg_noisy, glob_noisy] = v_snrseg(close_mic, noisy, 16e3, 'qbwp', 0.00625);
% seg_snr_improve = abs(seg_noisy-seg_enh)
% snr_improve = abs(glob_noisy-glob_enh)

figure();
v_sigalign(enh_irm, close_mic, 10, 'up', 16000);

figure();
snrseg(noisy, enhanced, 16e3, 'wp');

figure();
snrseg(noisy, enhanced2, 16e3, 'wp');

figure();
snrseg(enhanced, enhanced2, 16e3, 'wp');


fs = 16e3;
fftLength = 4096;
noverlap = 512;
win = hamming(fftLength,"periodic");

[audioIn,fs] = audioread("Laughter-16-8-mono-4secs.wav");
fs = 16e3;
[S,F,t] = stft(audioIn,fs, ...
               "Window",win, ...
               "OverlapLength",noverlap, ...
               "FFTLength",fftLength, ...
               "FrequencyRange","onesided");
PowerSpectrum = S.*conj(S);



numBands = 16;
normalization = "none";
[fb,cf] = designAuditoryFilterBank(32000, ...
            "FrequencyScale","erb", ...
            "FFTLength",fftLength, ...
            "NumBands",numBands, ...
            "FrequencyRange",[0 16000], ...
            "OneSided", true, ...
            "Normalization",normalization);
plot(F,fb.',LineWidth=2)
grid on
title("Mel Filter Bank")
xlabel("Frequency (Hz)")


bark_arr = zeros(4, size(f, 2));

[b1, c1] = v_frq2bark(f, 'lhg');
% [b11, c11] = v_frq2bark(f, 'lhgu');
[b2, c2] = v_frq2bark(f, 'LHg');
[b3, c3] = v_frq2bark(f, 'zg');
[b4, c4] = v_frq2bark(f, 'sg');

[mel1, mr1] = v_frq2mel(f);
v_frq2mel(f)

n_mels = 10;
n_fft = 1024;
fs = 16e3;
f_min = 0;
f_max = fs / 2;
scale = 'bg';

% f = 0:500:22e3;
% [erb, band] = v_frq2erb(f)
% ff = v_erb2frq(erb)
% 
% % x=v_filtbankm(p,n,fs,0,fs/2,'m');
% [x,cf]=v_filtbankm(n_mels, n_fft , fs, f_min, f_max, scale);
% 
% w = 'bg'
% 
% x=v_gammabank(n_mels, fs, w);
% 
% [b,a,fx,bx,gd,ph]=v_gammabank(0.35,fs,'g',[100 6000]);



% plot((0:floor(n_fft/2))*fs/n_fft,x')

bark_arr(1:4,:) = [b1; b2; b3; b4];
bark_arr_max = max(bark_arr, [], 1);
bark_arr_min = min(bark_arr, [], 1);
bark_max_dist = abs(bark_arr_max - bark_arr_min);
% v_frq2mel(f);

[f1, cf1] = v_bark2frq(b1, 'lhg');

figure();
subplot(3,1,1);
plot(f, bark_arr(1:4,:)', 'LineWidth',2)
grid on;
legend({'Traunmuller + corr', 'Traunmuller', 'Zwicker et al', 'Schroeder et al'});
xlim([0 22e3]);
xlabel('Frequency [Hz]');
ylabel('Bark Band [N]');
title('Bark-Band vs. Frequency');

subplot(3,1,2);
plot(f, bark_max_dist);
grid on;
xlim([0 22e3]);
xlabel('Frequency [Hz]');
ylabel('D_{max}');
title('Maximum Bark-Band Distance');

subplot(3,1,3);

bark_arr_max = max(bark_arr(1:3, :), [], 1);
bark_arr_min = min(bark_arr(1:3, :), [], 1);
bark_max_dist = abs(bark_arr_max - bark_arr_min);
plot(f, bark_max_dist);
grid on;
xlim([0 22e3]);
xlabel('Frequency [Hz]');
ylabel('D_{max}');
title('Maximum Bark-Band Distance, Exc. Bark4');

af = gcf;
exportgraphics(af,'bark_comparison2.jpg')


function [spectrums] = stft_voicebox(samples)
%using voicebox calculate multi-channel stft
frame_length = 400;
fft_length   = 512;
frame_shift  = 100;
hanning_wnd = v_windows(3,frame_length,'l')'; 
% hanning_wnd  = hanning(frame_length, 'periodic');
n_channel = size(samples, 2);

frames  = enframe(samples(:,1), hanning_wnd, frame_shift);
spectrums = zeros(fft_length/2+1, size(frames,1), n_channel);
for c = 1: n_channel
    frames  = enframe(samples(:,c), hanning_wnd, frame_shift);
    frames_padding = zeros(size(frames,1), fft_length);
    frames_padding(:, 1: frame_length) = frames;
    % rfft: T x F
    spectrums(:, :, c) = rfft(frames_padding, fft_length, 2)';
end
% spectrums = squeeze(spectrums);
end

function noise = estimate_noise_v2(clean, noisy)
%chime official method
wlen_sub=256; % STFT window length in samples
blen_sub=4000; % average block length in samples for speech subtraction (250 ms)
ntap_sub=12; % filter length in frames for speech subtraction (88 ms)
del=-3; % minimum delay (0 for a causal filter)
nsampl=length(clean);
% Compute the STFT (short window)
R=stft_multi(clean.',wlen_sub);
X=stft_multi(noisy.',wlen_sub);

% Estimate 88 ms impulse responses on 250 ms time blocks
A=estimate_ir(R,X,blen_sub,ntap_sub,del);

% Filter and subtract close-mic speech
Y=apply_ir(A,R,del);
y=istft_multi(Y,nsampl).';
noise=noisy-y;
end
% figure();
% subplot(2,1,1);
% plot(f, b1);
% hold on;
% plot(f, b2);
% grid on;
% plot(f, b3);
% plot(f, b4);
% legend;
% 
% subplot(2,1,2);
% plot(f, abs(abs(b2)-abs(b1)));


% figure();
% stem(c);
% grid on;
% legend;

% a = 10000;
% 
% x = 20:20e3;
% t1 = log(x);
% t2 = a*x.^(1/a)-a;
% % t3 = x.^2 / 2;
% 
% t4 = x.^(1/a);
% 
% figure();
% plot(t1);
% hold on;
% grid on;
% plot(t2);
% % plot(t3);
% % plot(t4);
% legend;