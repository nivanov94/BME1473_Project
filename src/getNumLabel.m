function class_num = getNumlabel(x)
% Function converts the label into nums
%   x  - Vector (1 x i) containing the original prompts i.e. '/iy/'
%   class_num   - Vector (1 x i) containing the numbered prompts i.e. '1'

class_word = x;
for i = 1:length(class_word)
    switch class_word{1,i}
        case '/iy/'
            class_num(1,i) = 1;
        case '/uw/'
            class_num(1,i) = 2;
        case '/piy/'
            class_num(1,i) = 3;
        case '/tiy/'
            class_num(1,i) = 4;
        case '/diy/'
            class_num(1,i) = 5;
        case '/m/'
            class_num(1,i) = 6;
        case '/n/'
            class_num(1,i) = 7;
        case 'pat'
            class_num(1,i) = 8;
        case 'pot'
            class_num(1,i) = 9;
        case 'knew'
            class_num(1,i) = 10;
        case 'gnaw'
            class_num(1,i) = 11;
        otherwise
            disp('unknown class exists');
            class_num(1,i) = 0;
    end
end
