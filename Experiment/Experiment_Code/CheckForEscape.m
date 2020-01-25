function CheckForEscape(kbID, keycode_esc, pahandle, res_FID)

if nargin == 3
    
    [keyisdown, firstpresstime] = KbQueueCheck(kbID);
    if firstpresstime(keycode_esc)
        KbQueueRelease; % Release keyboard internal queue
        PsychPortAudio('Close', pahandle); % Close the audio device
        sca % close screen
        Priority(0);
        ShowCursor;
        error('Escape key detected. Program exited!!!')
    end
end

if nargin == 4
    
    [keyisdown, firstpresstime] = KbQueueCheck(kbID);
    if firstpresstime(keycode_esc)
        fclose(res_FID);
        KbQueueRelease; % Release keyboard internal queue
        PsychPortAudio('Close', pahandle); % Close the audio device
        sca % close screen
        Priority(0);
        ShowCursor;
        error('Escape key detected. Program exited!!!')
    end
end

end