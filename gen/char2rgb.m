%%
%Name: char2rgb
%Created by: Pratik
%Created on: 27 July 2018
%
%Description: Converts a character color value into its numerical RGB equivalent
%Input: A lowercase character for any RGBCMYWK color
%Ouput: A 1x3 vector with the equivalent numerical RGB values
%
%Usage: vector = char2rgb('color');


%%
function rgbvec = char2rgb(color)
    switch color
        case 'r'
            rgbvec=[1 0 0];
        case 'g'
            rgbvec=[0 1 0];
        case 'b'
            rgbvec=[0 0 1];
        case 'c'
            rgbvec=[0 1 1];
        case 'm'
            rgbvec=[1 0 1];
        case 'y'
            rgbvec=[1 1 0];
        case 'w'
            rgbvec=[1 1 1];
        case 'k'
            rgbvec=[0 0 0];
        otherwise 
            return;
end