function[translated_state] = translate_states(input_state)
% Switch between notation for cell state:
% [0 0] = 1
% [1 0] = 2
% [0 1] = 3
% [1 1] = 4

if length(input_state) == 1
    if input_state == 1
        translated_state = [0 0];
    elseif input_state == 2
        translated_state = [1 0];
    elseif input_state == 3
        translated_state = [0 1];
    elseif input_state == 4
        translated_state = [1 1];
    else
        error('Inputted state not known.')
    end
else
    if input_state == [0 0]
        translated_state = 1;
    elseif input_state ==[1 0]
        translated_state = 2;
    elseif input_state == [0 1]
        translated_state = 3;
    elseif input_state == [1 1]
        translated_state = 4;
    else
        error('Inputted state not known.')
    end
end

end