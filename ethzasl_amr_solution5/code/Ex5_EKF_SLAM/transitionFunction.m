function [f, F_x, F_u] = transitionFunction(x,u,b)
% [f, F_x, F_u] = transitionFunction(x,u,b) predicts the state x at time t given
% the state at time t-1 and the input u at time t. b is the distance between
% the wheels f the differential-drive robot. F_x denotes the Jacobian
% of the state transition function with respect to the state evaluated at
% the state and input provided. F_u denotes the Jacobian of the state
% transition function with respect to the input evaluated at the state and
% input provided.
% State and input are defined according to the book pp. 337

if nargin < 3 || isempty(b)
    b = .1;
end

%STARTRM
f = x;
f(1) = f(1) + (u(1)+(u(2)))/2 * cos(x(3) + (u(2)-u(1))/(2*b) );
f(2) = f(2) + (u(1)+(u(2)))/2 * sin(x(3) + (u(2)-u(1))/(2*b) );
f(3) = f(3) + (u(2)-u(1))/(b);
%ENDRM
%STARTUNCOMMENT
%f(:) = TODO
%ENDUNCOMMENT

F_x = eye(size(x,1));
%STARTRM
F_x(1:3,1:3) = [...
    1, 0, -sin(x(3) - (u(1) - u(2))/(2*b))*(u(1)/2 + u(2)/2)
    0, 1,  cos(x(3) - (u(1) - u(2))/(2*b))*(u(1)/2 + u(2)/2)
    0, 0,  1
    ];
%ENDRM
%STARTUNCOMMENT
%F_x = TODO
%ENDUNCOMMENT

F_u = zeros(size(x,1),2);
%STARTRM
F_u(1:3,1:2) = [...
    cos(x(3) - (u(1) - u(2))/(2*b))/2 + (sin(x(3) - (u(1) - u(2))/(2*b))*(u(1)/2 + u(2)/2))/(2*b), cos(x(3) - (u(1) - u(2))/(2*b))/2 - (sin(x(3) - (u(1) - u(2))/(2*b))*(u(1)/2 + u(2)/2))/(2*b)
    sin(x(3) - (u(1) - u(2))/(2*b))/2 - (cos(x(3) - (u(1) - u(2))/(2*b))*(u(1)/2 + u(2)/2))/(2*b), sin(x(3) - (u(1) - u(2))/(2*b))/2 + (cos(x(3) - (u(1) - u(2))/(2*b))*(u(1)/2 + u(2)/2))/(2*b)
    -1/b,                                                                                      1/b
    ];
%ENDRM
%STARTUNCOMMENT
%F_u = TODO
%ENDUNCOMMENT

