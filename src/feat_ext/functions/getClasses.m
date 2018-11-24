function classes = getClasses( prompts )

% CLASSES ARE:
%       1: Vowel only (0) vs consonant (1)
%       2: non-nasal (0) vs nasal (1)
%       3: non-bilabial (0) vs bilabial (1)
%       4: non-iy (0) vs iy (1)
%       5: non-uw (0) vs uw (1)

class.iy  =     [ 0 0 0 1 0];
class.uw  =     [ 0 0 0 0 1];
class.m   =     [ 1 1 1 0 0];
class.n   =     [ 1 1 0 0 0];
class.piy =     [ 1 0 1 1 0];
class.tiy =     [ 1 0 0 1 0];
class.diy =     [ 1 0 0 1 0];
class.gnaw =    [ 1 1 0 0 0];
class.knew =    [ 1 1 0 0 0];
class.pat =     [ 1 0 1 0 0];
class.pot =     [ 1 0 1 0 0];

classes = [];

for i=1:length(prompts)
  prompt = regexprep( prompts{i}, '/', '');

  classes = [classes; class.(prompt)];
end