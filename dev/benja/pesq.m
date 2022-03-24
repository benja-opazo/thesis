function pesq_value = pesq(ref,deg) 
    %% FFMPEG 
    ffmpeg_output_name_ref = split(ref,'.wav');
    ffmpeg_output_name_deg = split(deg,'.wav');
 
    ffmpeg_command_ref = "ffmpeg -y -i " + ref + " -ar 16000 " + ffmpeg_output_name_ref{1} + '_lq.wav';
    ffmpeg_command_deg = "ffmpeg -y -i " + deg + " -ar 16000 " + ffmpeg_output_name_deg{1} + '_lq.wav';
    
    [~,~] = system(ffmpeg_command_ref);
    [~,~] = system(ffmpeg_command_deg);
    
    %% Python-PESQ
    
    %ref_file = '../outputs/tmp/fem_normal_aa_lowfs.wav';
    %deg_file = '../outputs/tmp/20210826-071120.wav';

    python_command = "python ../python-PESQ/pesq-python.py -r " + ffmpeg_output_name_ref{1} + '_lq.wav' + " -d " + ffmpeg_output_name_deg{1} + '_lq.wav';
    
    [~, pesq_value_str] = system(python_command);
    
    pesq_value = str2num(pesq_value_str);
    
    
    
end