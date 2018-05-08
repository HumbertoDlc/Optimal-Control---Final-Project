%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                         OPTIMAL CONTROL                         %%%%%
%%%%%                                                                 %%%%%
%%%%%                                                                 %%%%%
%%%%%                          Final Project                          %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Description of the states and inputs:
%%%%%
%%%%% 1: Ankle.
%%%%% 2: Knee.
%%%%% 3: Hip.
%%%%% 4: Shoulder.
%%%%% 5: Elbow.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% States and inputs:
syms q1 q2 q3 q4 q5 q1dot q2dot q3dot q4dot q5dot 
syms q1ddot q2ddot q3ddot q4ddot q5ddot u1 u2 u3 u4 u5 u6
U=[u1;u2;u3;u4;u5];                         % Input: Joint torques.
U(6)=u6;                                    % Input: Handle force.
x=[q1;q2;q3;q4;q5];                         % States: Joint position.
x(6:10)=[q1dot;q2dot;q3dot;q4dot;q5dot];    % States: Joint velocity.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Conditions:
q1=pi/4;    q1dot=0;    u1=0;
q2=pi/4;    q2dot=0;    u2=0;
q3=pi/4;    q3dot=0;    u3=0;
syms z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Body segment positions and lengths:       % (Units).
Mass        = 80;                               % Body mass (kg).
Xankle      = 78.8/1000;                        % Ankle x-position (mm).
Yankle      = 521.6/1000;                       % Ankle y-position (mm).
Xsprocket   = 1;                                % Sprocket x-position (mm).
Ysprocket   = 1;                                % Sprocket y-position (mm).
a           = 0.0286;          
b           = 786.6/1000;            
ShankLen    = 427.3/1000;                       % Shank length (mm).           
ThighLen    = 406.4/1000;                       % Thigh length (mm).
TrunkLen    = 428.7/1000;                       % Trunk length (mm).
UpparmLen   = 277.6/1000;                       % Upper arm length (mm).
ForearmLen  = 228.1/1000;                       % Forearm length (mm).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Body segment mass from Winter's book:     % (Units).
ShankMass   =  0.0465*2*Mass;                   % Shank mass (kg).
ThighMass   =  0.1*2*Mass;                      % Thigh mass (kg).
TrunkMass   =  0.497*Mass;                      % Trunk mass (kg).
UpparmMass  =  0.028*2*Mass;                    % Upper arm mass (kg).
ForearmMass =  0.016*2*Mass;                    % Forearm mass (kg).
ShankI      =  ShankMass*(0.302*ShankLen)^2;    % Shank inertia (kg.mm2).
ThighI      =  ThighMass*(0.323*ThighLen)^2;    % Thigh inertia (kg.mm2).
TrunkI      =  TrunkMass*(0.5*TrunkLen)^2;      % Trunk inertia (kg.mm2).
UpparmI     =  UpparmMass*(0.322*UpparmLen)^2;  % Upper arm inertia (kg.mm2).
ForearmI    =  ForearmMass*(0.303*ForearmLen)^2;% Forearm inertia (kg.mm2).   
ShankCM     =  (1-0.433)*ShankLen;              % Shank center-mass (mm).
ThighCM     =  (1-0.433)*ThighLen;              % Thigh center-mass (mm).
TrunkCM     =  (1-0.500)*TrunkLen;              % Trunk center-mass (mm).
UpparmCM    =  0.436*UpparmLen;                 % Upper center-mass (mm).
ForearmCM   =  0.430*ForearmLen;                % Forearm center-mass (mm).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Seat parameters and gravity:              % (Units).
Seat        =  0.1;                             % Seat mass (kg).	
Kseat       =  10^5;                            % Seat stiffness (kN/m).
Cseat       =  2*sqrt(Seat*Kseat);              % Seat damping (kNs/m).
g           =  9.80665;                         % gravity(N/kg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Dynamic equations based on the human physiology and states:
    z(1) = cos(q1);
    z(2) = sin(q1);
    z(5) = ShankLen*q1dot;
    z(6) = q1dot*z(5);
    z(7) = cos(q2);
    z(8) = sin(q2);
    z(9) = z(1)*z(7) - z(2)*z(8);
    z(10) = -z(1)*z(8) - z(2)*z(7);
    z(11) = z(1)*z(8) + z(2)*z(7);
    z(12) = ShankLen*z(8);
    z(13) = ShankLen*z(7);
    z(14) = ThighCM + z(13);
    z(15) = z(7)*z(6);
    z(16) = z(8)*z(6);
    z(17) = ThighCM*(q1dot+q2dot);
    z(18) = q1dot + q2dot;
    z(19) = -z(15) - z(17)*z(18);
    z(20) = ThighLen + z(13);
    z(21) = ThighLen*(q1dot+q2dot);
    z(22) = -z(15) - z(18)*z(21);
    z(23) = cos(q3);
    z(24) = sin(q3);
    z(25) = z(12)*z(23) + z(20)*z(24);
    z(26) = ThighLen*z(24);
    z(27) = ThighLen*z(23);
    z(28) = z(20)*z(23) - z(12)*z(24);
    z(29) = z(23)^2 + z(24)^2;
    z(30) = TrunkCM*(1+z(29));
    z(31) = z(27) - TrunkCM;
    z(32) = z(28) - TrunkCM;
    z(33) = z(23)*z(22) + z(24)*z(16);
    z(34) = z(23)*z(16) - z(24)*z(22);
    z(36) = q1dot + q2dot + q3dot;
    z(37) = z(29)*q3dot;
    z(38) = z(33) + (TrunkCM*q1dot+TrunkCM*q2dot+z(30)*q3dot)*(z(36)+z(37));
    z(39) = TrunkLen*(1+z(29));
    z(40) = z(27) - TrunkLen;
    z(41) = z(28) - TrunkLen;
    z(42) = z(33) + (TrunkLen*q1dot+TrunkLen*q2dot+z(39)*q3dot)*(z(36)+z(37));
    z(43) = cos(q4);
    z(44) = sin(q4);
    z(45) = z(9)*z(23) + z(10)*z(24);
    z(46) = z(10)*z(23) - z(9)*z(24);
    z(47) = z(9)*z(24) + z(11)*z(23);
    z(48) = z(9)*z(23) - z(11)*z(24);
    z(49) = z(43)*z(45) + z(44)*z(46);
    z(50) = z(43)*z(46) - z(44)*z(45);
    z(51) = z(43)*z(47) + z(44)*z(48);
    z(52) = z(43)*z(48) - z(44)*z(47);
    z(53) = z(25)*z(43) + z(41)*z(44);
    z(54) = z(26)*z(43) + z(40)*z(44);
    z(55) = z(39)*z(44);
    z(56) = z(41)*z(43) - z(25)*z(44);
    z(57) = z(40)*z(43) - z(26)*z(44);
    z(58) = z(39)*z(43);
    z(59) = UpparmCM + z(56);
    z(60) = UpparmCM + z(57);
    z(61) = UpparmCM - z(58);
    z(62) = z(43)*z(42) + z(44)*z(34);
    z(63) = z(43)*z(34) - z(44)*z(42);
    z(64) = UpparmCM*(q1dot+q2dot+q3dot+q4dot);
    z(65) = q1dot + q2dot + q3dot + q4dot;
    z(66) = z(62) - z(64)*z(65);
    z(67) = UpparmLen + z(56);
    z(68) = UpparmLen + z(57);
    z(69) = UpparmLen - z(58);
    z(70) = UpparmLen*(q1dot+q2dot+q3dot+q4dot);
    z(71) = z(62) - z(65)*z(70);
    z(72) = cos(q5);
    z(73) = sin(q5);
    z(74) = z(49)*z(72) + z(50)*z(73);
    z(75) = z(50)*z(72) - z(49)*z(73);
    z(76) = z(51)*z(72) + z(52)*z(73);
    z(77) = z(52)*z(72) - z(51)*z(73);
    z(78) = z(53)*z(72) + z(67)*z(73);
    z(79) = z(54)*z(72) + z(68)*z(73);
    z(80) = z(69)*z(73) - z(55)*z(72);
    z(81) = UpparmLen*z(73);
    z(82) = UpparmLen*z(72);
    z(83) = z(67)*z(72) - z(53)*z(73);
    z(84) = z(68)*z(72) - z(54)*z(73);
    z(85) = z(55)*z(73) + z(69)*z(72);
    z(86) = ForearmCM + z(82);
    z(87) = ForearmCM + z(83);
    z(88) = ForearmCM + z(84);
    z(89) = ForearmCM + z(85);
    z(90) = z(72)*z(71) + z(73)*z(63);
    z(91) = z(72)*z(63) - z(73)*z(71);
    z(92) = ForearmCM*(q1dot+q2dot+q3dot+q4dot+q5dot);
    z(93) = q1dot + q2dot + q3dot + q4dot + q5dot;
    z(94) = z(90) - z(92)*z(93);
    z(95) = ForearmMass*g;
    z(96) = g*ShankMass;
    z(97) = g*ThighMass;
    z(98) = g*TrunkMass;
    z(99) = g*UpparmMass;
    d     = -(Yankle+ShankLen*z(2)+ThighLen*z(11)-b-a*Xankle...
            -a*ShankLen*z(1)-a*ThighLen*z(9))/(1+a^2)^0.5;
    z(100) = a*ShankLen;
    z(101) = a*ThighLen;
    z(102) = 1 + a^2;
    z(103) = (ShankLen*z(1)+ThighLen*z(9)+z(100)*z(2)-...
             z(101)*z(10))/z(102)^0.5;
    z(104) = ThighLen*(z(9)-a*z(10))/z(102)^0.5;
    z(105) = a*(Kseat*d-Cseat*(z(103)*q1dot+z(104)*q2dot))/z(102)^0.5;
    z(106) = (Kseat*d-Cseat*(z(103)*q1dot+z(104)*q2dot))/z(102)^0.5;
    z(107) = ForearmLen + z(82);
    z(108) = ForearmLen + z(83);
    z(109) = ForearmLen + z(84);
    z(110) = ForearmLen + z(85);
    z(111) = Xsprocket - Xankle;
    z(112) = Ysprocket - Yankle;
    z(113) = z(1)*z(74) + z(2)*z(76);
    z(115) = z(1)*z(75) + z(2)*z(77);
    z(117) = z(9)*z(74) + z(11)*z(76);
    z(119) = z(9)*z(75) + z(11)*z(77);
    z(121) = z(43)*z(72) - z(44)*z(73);
    z(123) = -z(43)*z(73) - z(44)*z(72);
    z(124) = z(1)*z(45) + z(2)*z(47);
    z(125) = z(1)*z(46) + z(2)*z(48);
    z(128) = z(1)*z(49) + z(2)*z(51);
    z(129) = z(1)*z(50) + z(2)*z(52);
    z(130) = z(1)*z(51) - z(2)*z(49);
    z(131) = z(1)*z(52) - z(2)*z(50);
    z(132) = z(7)*z(128) + z(8)*z(130);
    z(133) = z(7)*z(129) + z(8)*z(131);
    L      = (ForearmLen^2+ShankLen^2+ThighLen^2+TrunkLen^2+UpparmLen^2+...
             z(111)^2+z(112)^2+2*ForearmLen*ShankLen*z(113)+2*...
             ForearmLen*ThighLen*z(117)+2*ForearmLen*UpparmLen*z(72)+...
             2*ShankLen*ThighLen*z(7)+2*ShankLen*UpparmLen*z(128)+2*...
             ThighLen*UpparmLen*z(132)+2*TrunkLen*z(111)*z(45)+2*...
             TrunkLen*z(112)*z(47)-2*ForearmLen*TrunkLen*z(121)-2*...
             ForearmLen*z(111)*z(74)-2*ForearmLen*z(112)*z(76)-2*...
             ShankLen*TrunkLen*z(124)-2*ShankLen*z(111)*z(1)-2*...
             ShankLen*z(112)*z(2)-2*ThighLen*TrunkLen*z(23)-2*ThighLen*...
             z(111)*z(9)-2*ThighLen*z(112)*z(11)-2*TrunkLen*UpparmLen*...
             z(43)-2*UpparmLen*z(111)*z(49)-2*UpparmLen*z(112)*z(51))^0.5;
    z(136) = z(2)*z(8) - z(1)*z(7);
    z(137) = z(10)*z(23) + z(24)*z(136);
    z(138) = z(23)*z(136) - z(10)*z(24);
    z(139) = -z(9)*z(23) - z(10)*z(24);
    z(140) = z(43)*z(46) + z(44)*z(139);
    z(141) = z(43)*z(137) + z(44)*z(138);
    z(142) = z(43)*z(138) - z(44)*z(137);
    z(143) = z(43)*z(139) - z(44)*z(46);
    z(144) = -z(43)*z(45) - z(44)*z(46);
    z(145) = z(50)*z(72) + z(73)*z(144);
    z(146) = z(72)*z(140) + z(73)*z(143);
    z(147) = z(72)*z(141) + z(73)*z(142);
    z(148) = -z(9)*z(24) - z(11)*z(23);
    z(149) = z(43)*z(48) + z(44)*z(148);
    z(150) = z(43)*z(148) - z(44)*z(48);
    z(151) = -z(43)*z(47) - z(44)*z(48);
    z(152) = z(52)*z(72) + z(73)*z(151);
    z(153) = z(72)*z(149) + z(73)*z(150);
    z(154) = z(1)*(z(76)+z(147));
    z(155) = z(1)*z(145) + z(2)*z(152);
    z(156) = z(1)*z(146) + z(2)*z(153);
    z(157) = z(1)*z(147) + z(2)*z(74);
    z(158) = z(9)*z(76) + z(9)*z(147) + z(10)*z(74) + z(11)*z(74);
    z(159) = z(9)*z(145) + z(11)*z(152);
    z(160) = z(9)*z(146) + z(11)*z(153);
    z(161) = z(1)*(z(51)+z(141));
    z(162) = z(1)*z(140) + z(2)*z(149);
    z(163) = z(1)*z(141) + z(2)*z(49);
    z(164) = z(2)*(z(51)+z(141));
    z(165) = z(1)*z(49) - z(2)*z(141);
    z(166) = z(1)*z(149) - z(2)*z(140);
    z(167) = z(7)*z(130) + z(7)*z(163) + z(8)*z(165) - z(8)*z(128);
    z(168) = z(7)*z(161) - z(8)*z(164);
    z(169) = z(7)*z(162) + z(8)*z(166);
    z(170) = z(1)*(z(47)+z(137));
    z(171) = z(1)*z(137) + z(2)*z(45);
    z(172) = ForearmLen*ShankLen;
    z(173) = ForearmLen*ThighLen;
    z(174) = ShankLen*UpparmLen;
    z(175) = ShankLen*z(111);
    z(176) = ThighLen*UpparmLen;
    z(177) = TrunkLen*z(111);
    z(178) = TrunkLen*z(112);
    z(179) = ForearmLen*z(111);
    z(180) = ForearmLen*z(112);
    z(181) = ShankLen*TrunkLen;
    z(182) = ShankLen*z(112);
    z(183) = ThighLen*z(111);
    z(184) = ThighLen*z(112);
    z(185) = UpparmLen*z(111);
    z(186) = UpparmLen*z(112);
    z(187) = ForearmLen^2 + ShankLen^2 + ThighLen^2 + TrunkLen^2 + ...
             UpparmLen^2 + z(111)^2 + z(112)^2;
    z(188) = ForearmLen*UpparmLen;
    z(189) = ShankLen*ThighLen;
    z(190) = ForearmLen*TrunkLen;
    z(191) = ThighLen*TrunkLen;
    z(192) = TrunkLen*UpparmLen;
    z(193) = (z(172)*z(154)+z(173)*z(158)+z(174)*z(161)+z(175)*...
             z(2)+z(176)*z(168)+z(177)*z(137)+z(178)*z(45)-z(179)*...
             z(147)-z(180)*z(74)-z(181)*z(170)-z(182)*z(1)-z(183)*...
             z(10)-z(184)*z(9)-z(185)*z(141)-z(186)*z(49))/(z(187)+2*...
             z(172)*z(113)+2*z(173)*z(117)+2*z(174)*z(128)+2*z(176)*...
             z(132)+2*z(177)*z(45)+2*z(178)*z(47)+2*z(188)*z(72)+2*...
             z(189)*z(7)-2*z(175)*z(1)-2*z(179)*z(74)-2*z(180)*z(76)-2*...
             z(181)*z(124)-2*z(182)*z(2)-2*z(183)*z(9)-2*z(184)*z(11)-2*...
             z(185)*z(49)-2*z(186)*z(51)-2*z(190)*z(121)-2*z(191)*z(23)-...
             2*z(192)*z(43))^0.5;
    z(194) = (z(179)*z(146)+z(180)*z(153)+z(181)*z(125)+z(185)*z(140)+...
             z(186)*z(149)-z(172)*z(156)-z(173)*z(160)-z(174)*z(162)-...
             z(176)*z(169)-z(177)*z(46)-z(178)*z(48)-z(191)*z(24))/...
             (z(187)+2*z(172)*z(113)+2*z(173)*z(117)+2*z(174)*z(128)+...
             2*z(176)*z(132)+2*z(177)*z(45)+2*z(178)*z(47)+2*z(188)*...
             z(72)+2*z(189)*z(7)-2*z(175)*z(1)-2*z(179)*z(74)-2*z(180)*...
             z(76)-2*z(181)*z(124)-2*z(182)*z(2)-2*z(183)*z(9)-2*z(184)*...
             z(11)-2*z(185)*z(49)-2*z(186)*z(51)-2*z(190)*z(121)-2*...
             z(191)*z(23)-2*z(192)*z(43))^0.5;
    z(195) = (z(179)*z(145)+z(180)*z(152)+z(185)*z(50)+z(186)*z(52)+...
             z(190)*z(123)-z(172)*z(155)-z(173)*z(159)-z(174)*z(129)-...
             z(176)*z(133)-z(192)*z(44))/(z(187)+2*z(172)*z(113)+2*...
             z(173)*z(117)+2*z(174)*z(128)+2*z(176)*z(132)+2*z(177)*...
             z(45)+2*z(178)*z(47)+2*z(188)*z(72)+2*z(189)*z(7)-2*...
             z(175)*z(1)-2*z(179)*z(74)-2*z(180)*z(76)-2*z(181)*...
             z(124)-2*z(182)*z(2)-2*z(183)*z(9)-2*z(184)*z(11)-2*...
             z(185)*z(49)-2*z(186)*z(51)-2*z(190)*z(121)-2*z(191)*...
             z(23)-2*z(192)*z(43))^0.5;
    z(196) = (z(172)*z(157)+z(173)*z(158)+z(174)*z(163)+z(176)*...
             z(167)+z(177)*z(137)+z(178)*z(45)-z(179)*z(147)-z(180)*...
             z(74)-z(181)*z(171)-z(183)*z(10)-z(184)*z(9)-z(185)*...
             z(141)-z(186)*z(49)-z(189)*z(8))/(z(187)+2*z(172)*...
             z(113)+2*z(173)*z(117)+2*z(174)*z(128)+2*z(176)*z(132)+...
             2*z(177)*z(45)+2*z(178)*z(47)+2*z(188)*z(72)+2*z(189)*...
             z(7)-2*z(175)*z(1)-2*z(179)*z(74)-2*z(180)*z(76)-2*...
             z(181)*z(124)-2*z(182)*z(2)-2*z(183)*z(9)-2*z(184)*...
             z(11)-2*z(185)*z(49)-2*z(186)*z(51)-2*z(190)*z(121)-2*...
             z(191)*z(23)-2*z(192)*z(43))^0.5;
    z(197) = ForearmLen*(ShankLen*z(115)+ThighLen*z(119)-TrunkLen*...
             z(123)-UpparmLen*z(73)-z(111)*z(75)-z(112)*z(77))/(z(187)+...
             2*z(172)*z(113)+2*z(173)*z(117)+2*z(174)*z(128)+...
             2*z(176)*z(132)+2*z(177)*z(45)+2*z(178)*z(47)+2*z(188)*...
             z(72)+2*z(189)*z(7)-2*z(175)*z(1)-2*z(179)*z(74)-2*z(180)*...
             z(76)-2*z(181)*z(124)-2*z(182)*z(2)-2*z(183)*z(9)-2*z(184)*...
             z(11)-2*z(185)*z(49)-2*z(186)*z(51)-2*z(190)*z(121)-2*...
             z(191)*z(23)-2*z(192)*z(43))^0.5;
    z(198) = z(187) + 2*z(172)*z(113) + 2*z(173)*z(117) + 2*z(174)*...
             z(128) + 2*z(176)*z(132) + 2*z(177)*z(45) + 2*z(178)*z(47)+... 
             2*z(188)*z(72) + 2*z(189)*z(7) - 2*z(175)*z(1) - 2*z(179)*...
             z(74) - 2*z(180)*z(76) - 2*z(181)*z(124) - 2*z(182)*z(2) - ...
             2*z(183)*z(9) - 2*z(184)*z(11) - 2*z(185)*z(49) - 2*z(186)*...
             z(51) - 2*z(190)*z(121) - 2*z(191)*z(23) - 2*z(192)*z(43);
    z(199) = z(198)^0.5;
    z(200) = ForearmLen/z(199);
    z(201) = z(111)/z(199);
    z(202) = z(112)/z(199);
    z(203) = ShankLen/z(199);
    z(204) = ThighLen/z(199);
    z(205) = TrunkLen/z(199);
    z(206) = UpparmLen/z(199);
    z(207) = U(6)*z(200);
    z(208) = U(6)*z(201);
    z(209) = U(6)*z(202);
    z(210) = U(6)*z(203);
    z(211) = U(6)*z(204);
    z(212) = U(6)*z(205);
    z(213) = U(6)*z(206);
    z(214) = ShankCM*z(96);
    z(215) = z(73)*z(108)*z(213) + z(74)*z(78)*z(208) + z(75)*z(108)*...
             z(208) + z(76)*z(78)*z(209) + z(77)*z(108)*z(209) + z(78)*...
             z(121)*z(212) + z(108)*z(123)*z(212) + z(9)*z(20)*z(106) +...
             z(11)*z(12)*z(106) - z(214)*z(1) - z(78)*z(207) - z(72)*...
             z(78)*z(213) - z(78)*z(113)*z(210) - z(78)*z(117)*z(211) -...
             z(108)*z(115)*z(210) - z(108)*z(119)*z(211) - z(95)*(z(76)*...
             z(78)+z(77)*z(87)) - z(97)*(z(9)*z(14)+z(11)*z(12)) - ...
             z(98)*(z(25)*z(47)+z(32)*z(48)) - z(99)*(z(51)*z(53)+z(52)*...
             z(59)) - z(9)*z(12)*z(105) - z(10)*z(20)*z(105);
    z(216) = ThighCM*z(97);
    z(217) = z(73)*z(109)*z(213) + z(74)*z(79)*z(208) + z(75)*z(109)*...
             z(208) + z(76)*z(79)*z(209) + z(77)*z(109)*z(209) + z(79)*...
             z(121)*z(212) + z(109)*z(123)*z(212) + ThighLen*(z(9)*...
             z(106)-z(10)*z(105)) - z(216)*z(9) - z(79)*z(207) - z(72)*...
             z(79)*z(213) - z(79)*z(113)*z(210) - z(79)*z(117)*z(211) -...
             z(109)*z(115)*z(210) - z(109)*z(119)*z(211) - z(95)*(z(76)*...
             z(79)+z(77)*z(88)) - z(98)*(z(26)*z(47)+z(31)*z(48)) - ...
             z(99)*(z(51)*z(54)+z(52)*z(60));
    z(218) = z(98)*z(30)*z(48) + z(73)*z(110)*z(213) + z(74)*z(80)*...
             z(208) + z(75)*z(110)*z(208) + z(76)*z(80)*z(209) + z(77)*...
             z(110)*z(209) + z(80)*z(121)*z(212) + z(110)*z(123)*z(212)+...
             z(99)*(z(51)*z(55)-z(52)*z(61)) - z(80)*z(207) - z(72)*...
             z(80)*z(213) - z(80)*z(113)*z(210) - z(80)*z(117)*z(211) - ...
             z(110)*z(115)*z(210) - z(110)*z(119)*z(211) - z(95)*(z(76)*...
             z(80)+z(77)*z(89));
    z(219) = UpparmCM*z(99);
    z(220) = z(73)*z(107)*z(213) + z(74)*z(81)*z(208) + z(75)*z(107)*...
             z(208) + z(76)*z(81)*z(209) + z(77)*z(107)*z(209) + z(81)*...
             z(121)*z(212) + z(107)*z(123)*z(212) - z(219)*z(52) - ...
             z(81)*z(207) - z(72)*z(81)*z(213) - z(81)*z(113)*z(210) - ...
             z(81)*z(117)*z(211) - z(107)*z(115)*z(210) - z(107)*...
             z(119)*z(211) - z(95)*(z(76)*z(81)+z(77)*z(86));
    z(221) = ForearmCM*z(95);
    z(222) = -z(221)*z(77) - ForearmLen*(z(115)*z(210)+z(119)*z(211)-...
             z(73)*z(213)-z(75)*z(208)-z(77)*z(209)-z(123)*z(212));
    z(223) = z(45)*z(48) + z(46)*z(137);
    z(224) = z(46)^2 + z(48)^2;
    z(226) = -z(1)*z(8)*q1dot - z(1)*z(8)*q2dot - z(2)*z(7)*q1dot -...
             z(2)*z(7)*q2dot;
    z(227) = z(2)*z(8)*q1dot + z(2)*z(8)*q2dot - z(1)*z(7)*q1dot - ...
             z(1)*z(7)*q2dot;
    z(228) = z(10)*z(23)*q3dot+z(23)*z(226)+z(24)*z(227)-z(9)*z(24)*q3dot;
    z(229) = z(1)*z(7)*q1dot + z(1)*z(7)*q2dot - z(2)*z(8)*q1dot -...
             z(2)*z(8)*q2dot;
    z(230) = z(23)*z(226)-z(9)*z(24)*q3dot-z(11)*z(23)*q3dot-z(24)*z(229);
    z(231) = z(23)*z(227) - z(9)*z(23)*q3dot - z(10)*z(24)*q3dot - z(24)*z(226);
    z(232) = z(1)*z(8)*q1dot + z(1)*z(8)*q2dot + z(2)*z(7)*q1dot +...
             z(2)*z(7)*q2dot;
    z(233) = z(23)*z(136)*q3dot + z(23)*z(227) + z(24)*z(232) -...
             z(10)*z(24)*q3dot;
    z(234) = z(45)*q1dot*z(230) + z(45)*q2dot*z(230) + z(46)*q1dot*...
             z(233) + z(46)*q2dot*z(233) + z(48)*q1dot*z(228) + ...
             z(48)*q2dot*z(228) + z(137)*q1dot*z(231) + z(137)*q2dot*...
             z(231) + 2*z(46)*q3dot*z(231) + 2*z(48)*q3dot*z(230);
    z(235) = TrunkI*z(223);
    z(236) = TrunkI*z(224);
    z(237) = TrunkI*z(234);
    z(238) = ForearmI + ShankI + ThighI + UpparmI + ShankMass*ShankCM^2;
    z(239) = z(238) + z(223)*z(235) + ForearmMass*(z(78)^2+z(87)^2) +...
             ThighMass*(z(12)^2+z(14)^2) + TrunkMass*(z(25)^2+z(32)^2) +...
             UpparmMass*(z(53)^2+z(59)^2);
    z(240) = ForearmI + ThighI + UpparmI;
    z(241) = ThighCM*ThighMass;
    z(242) = z(240) + z(241)*z(14) + z(223)*z(235) + ForearmMass*(z(78)*...
             z(79)+z(87)*z(88)) + TrunkMass*(z(25)*z(26)+z(31)*z(32)) +...
             UpparmMass*(z(53)*z(54)+z(59)*z(60));
    z(243) = z(240) + z(223)*z(236) + ForearmMass*(z(78)*z(80)+z(87)*...
             z(89)) - TrunkMass*z(30)*z(32) - UpparmMass*(z(53)*z(55)-...
             z(59)*z(61));
    z(244) = ForearmI + UpparmI;
    z(245) = UpparmCM*UpparmMass;
    z(246) = z(244) + z(245)*z(59) + ForearmMass*(z(78)*z(81)+z(86)*z(87));
    z(247) = ForearmCM*ForearmMass;
    z(248) = ForearmI + z(247)*z(87);
    z(249) = z(223)*z(237) + ForearmMass*(z(78)*z(94)+z(87)*z(91)) +...
             ThighMass*(z(12)*z(19)+z(14)*z(16)) + TrunkMass*(z(25)*...
             z(38)+z(32)*z(34)) + UpparmMass*(z(53)*z(66)+z(59)*z(63));
    z(250) = ForearmI + ThighI + UpparmI + ThighMass*ThighCM^2;
    z(251) = z(250) + z(223)*z(235) + ForearmMass*(z(79)^2+z(88)^2) +...
             TrunkMass*(z(26)^2+z(31)^2) + UpparmMass*(z(54)^2+z(60)^2);
    z(252) = z(240) + z(223)*z(236) + ForearmMass*(z(79)*z(80)+z(88)*...
             z(89)) - TrunkMass*z(30)*z(31) - UpparmMass*(z(54)*z(55)-...
             z(60)*z(61));
    z(253) = z(244) + z(245)*z(60) + ForearmMass*(z(79)*z(81)+z(86)*z(88));
    z(254) = ForearmI + z(247)*z(88);
    z(255) = z(241)*z(16) + z(223)*z(237) + ForearmMass*(z(79)*z(94)+...
             z(88)*z(91)) + TrunkMass*(z(26)*z(38)+z(31)*z(34)) + ...
             UpparmMass*(z(54)*z(66)+z(60)*z(63));
    z(256) = z(240) + z(224)*z(235) + ForearmMass*(z(78)*z(80)+...
             z(87)*z(89)) - TrunkMass*z(30)*z(32) - UpparmMass*(z(53)*...
             z(55)-z(59)*z(61));
    z(257) = z(240) + z(224)*z(235) + ForearmMass*(z(79)*z(80)+z(88)*...
             z(89)) - TrunkMass*z(30)*z(31) - UpparmMass*(z(54)*z(55)-...
             z(60)*z(61));
    z(258) = z(240) + z(224)*z(236) + TrunkMass*z(30)^2 + ForearmMass*...
             (z(80)^2+z(89)^2) + UpparmMass*(z(55)^2+z(61)^2);
    z(259) = z(244) + z(245)*z(61) + ForearmMass*(z(80)*z(81)+z(86)*z(89));
    z(260) = ForearmI + z(247)*z(89);
    z(261) = z(224)*z(237) + ForearmMass*(z(80)*z(94)+z(89)*z(91)) -...
             TrunkMass*z(30)*z(34) - UpparmMass*(z(55)*z(66)-z(61)*z(63));
    z(262) = ForearmI + UpparmI + UpparmMass*UpparmCM^2;
    z(263) = z(262) + ForearmMass*(z(81)^2+z(86)^2);
    z(264) = ForearmI + z(247)*z(86);
    z(265) = z(245)*z(63) + ForearmMass*(z(81)*z(94)+z(86)*z(91));
    z(266) = ForearmI + ForearmMass*ForearmCM^2;
    z(267) = z(247)*z(91);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Dynamics parameters:
Inertia    = [      28.37,  0,      0,      0,      0;
                    0,      6.46,   0,      0,      0;
                    0,      0,      13.03,  0,      0;
                    0,      0,      0,      0.48,   0;
                    0,      0,      0,      0,      0.04];
Theta_1       = [   z(249);
                    z(255);
                    z(261);
                    z(265);
                    z(267)];

Theta_2       = [   z(215);
                    z(217);
                    z(218);
                    z(220);
                    z(222)];

Theta       =   Theta_1-Theta_2;    
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% State equations:
x_dot=x(6:10);
x_dot(6:10)=Inertia\(U(1:5)-Theta);
x_dot=vpa(simplify(x_dot),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cost function parameters: J= h + int(g)dt.
h=0;
g=5*(150-U(4))^2+5*(150-U(5))^2+U(6)^2+50*q4dot^2+50*q5dot^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Hamiltonian:
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10
H=g+p1*(x_dot(1))+p2*(x_dot(2))+p3*(x_dot(3))+p4*(x_dot(4))+p5*(x_dot(5))+...
    p6*(x_dot(6))+p7*(x_dot(7))+p8*(x_dot(8))+p9*(x_dot(9))+p10*(x_dot(10));
p1_dot=0;
p2_dot=0;
p3_dot=0;
p4_dot=vpa(simplify(-diff(H,q4)),2);
p5_dot=vpa(simplify(-diff(H,q5)),2);
p6_dot=0;
p7_dot=0;
p8_dot=0;
p9_dot=vpa(simplify(-diff(H,q4dot)),2);
p10_dot=vpa(simplify(-diff(H,q5dot)),2);
u_star_1=0;
u_star_2=0;
u_star_3=0;
u_star_4=vpa(simplify(diff(H,u4)),2);
u_star_5=vpa(simplify(diff(H,u5)),2);
u_star_6=vpa(simplify(diff(H,u6)),2);
