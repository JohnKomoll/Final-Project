% Project: Asthma and Allergies
%  Genes: ORMDL3, GSDML, IL-33, TSLP, RHO (Control)

% First, examine the sequences of the genes with BLAST

% Give accession numbers
acc_ormdl3 = 'NM_139280';
acc_gsmdl = 'BX538068';
acc_il33 = 'NM_033439';
acc_tslp = 'AF338732';

% Get genbank data
ormdl3_dat = getgenbank(acc_ormdl3);
gsmdl_dat = getgenbank(acc_gsmdl);
il33_dat = getgenbank(acc_il33);
tslp_dat = getgenbank(acc_tslp);

% Perform Smith-Waterman aligning algorithm on each pair of sequences
swalign(ormdl3_dat.Sequence, gsmdl_dat.Sequence, 'Showscore', true);
title('ORMDL3 and GSDML Sequence Comparison')
xlabel('ORMDL3 Sequence')
ylabel('GSDML Sequence')
swalign(ormdl3_dat.Sequence, il33_dat.Sequence, 'Showscore', true);
title('ORMDL3 and IL-33 Sequence Comparison')
xlabel('ORMDL3 Sequence')
ylabel('IL-33 Sequence')
swalign(ormdl3_dat.Sequence, tslp_dat.Sequence, 'Showscore', true);
title('ORMDL3 and TSLP Sequence Comparison')
xlabel('ORMDL3 Sequence')
ylabel('TSLP Sequence')
swalign(gsmdl_dat.Sequence, il33_dat.Sequence, 'Showscore', true);
title('GSDML and IL-33 Sequence Comparison')
xlabel('GSDML Sequence')
ylabel('IL-33 Sequence')
swalign(gsmdl_dat.Sequence, tslp_dat.Sequence, 'Showscore', true);
title('GSDML and TSLP Sequence Comparison')
xlabel('GSDML Sequence')
ylabel('TSLP Sequence')
swalign(il33_dat.Sequence, tslp_dat.Sequence, 'Showscore', true);
title('IL-33 and TSLP Sequence Comparison')
xlabel('IL-33 Sequence')
ylabel('TSLP Sequence')

%%
% Access BLAST data for each gene sequence
[reqID, reqtime] = blastncbi(ormdl3_dat.Sequence, 'blastn');
ormdl3_blast = getblast(reqID, 'WaitTime', reqtime);

[reqID, reqtime] = blastncbi(gsmdl_dat.Sequence, 'blastn');
gsmdl_blast = getblast(reqID, 'WaitTime', reqtime);

[reqID, reqtime] = blastncbi(il33_dat.Sequence, 'blastn');
il33_blast = getblast(reqID, 'WaitTime', reqtime);

[reqID, reqtime] = blastncbi(tslp_dat.Sequence, 'blastn');
tslp_blast = getblast(reqID, 'WaitTime', reqtime);

% All BLAST queries between each set of two sequences returned no
% significant alignment. This is probably due to the very complicated
% pathways that cause allergic and asthma-related reactions.

%%
% Check if they alter cell division at all?
% Try to fit 
% addpath 'C:\Users\John\Documents\BIOC 470 project\bfmatlab'

% Import videos
il33_vid = VideoReader('IL33/130--11--10--(14,19)--131654--C9orf26-gfp.mp4');
il33_vid2 = VideoReader('IL33/130--11--10--(14,19)--131654--C9orf26-gfp (1).mp4');
il33_vid3 = VideoReader('IL33/130--11--10--(14,19)--131654--C9orf26-gfp (2).mp4');
il33_vid4 = VideoReader('IL33/144--12--12--(16,19)--131655--C9orf26-gfp.mp4');
il33_vid5 = VideoReader('IL33/144--12--12--(16,19)--131655--C9orf26-gfp (1).mp4');
ormdl3_vid = VideoReader('ORDML3/195--17--03--(1,21)--129213--ORMDL3-gfp.mp4');
ormdl3_vid2 = VideoReader('ORDML3/195--17--03--(1,21)--129213--ORMDL3-gfp (1).mp4');
ormdl3_vid3 = VideoReader('ORDML3/195--17--03--(1,21)--129213--ORMDL3-gfp (2).mp4');
ormdl3_vid4 = VideoReader('ORDML3/217--19--01--(3,21)--129215--ORMDL3-gfp.mp4');
ormdl3_vid5 = VideoReader('ORDML3/217--19--01--(3,21)--129215--ORMDL3-gfp (1).mp4');
gsdml_vid = VideoReader('GSMDL/041--04--05--(10,9)--132569--GSDML-gfp.mp4');
gsdml_vid2 = VideoReader('GSMDL/041--04--05--(10,9)--132569--GSDML-gfp (1).mp4');
gsdml_vid3 = VideoReader('GSMDL/041--04--05--(10,9)--132569--GSDML-gfp (2).mp4');
gsdml_vid4 = VideoReader('GSMDL/011--01--11--(12,9)--132570--GSDML-gfp.mp4');
gsdml_vid5 = VideoReader('GSMDL/011--01--11--(12,9)--132570--GSDML-gfp (1).mp4');
tslp_vid = VideoReader('TSLP/373--32--01--(5,8)--127709--TSLP-gfp.mp4');
tslp_vid2 = VideoReader('TSLP/373--32--01--(5,8)--127709--TSLP-gfp (1).mp4');
tslp_vid3 = VideoReader('TSLP/373--32--01--(5,8)--127709--TSLP-gfp (2).mp4');
tslp_vid4 = VideoReader('TSLP/343--29--07--(7,8)--127710--TSLP-gfp.mp4');
tslp_vid5 = VideoReader('TSLP/343--29--07--(7,8)--127710--TSLP-gfp (1).mp4');
control_vid = VideoReader('RHO/107--09--11--(12,11)--2171--RHO-gfp.mp4');
control_vid2 = VideoReader('RHO/107--09--11--(12,11)--2171--RHO-gfp (1).mp4');
control_vid3 = VideoReader('RHO/229--20--01--(5,5)--2080--RHO-gfp.mp4');
control_vid4 = VideoReader('RHO/296--25--08--(7,15)--42859--RHO-gfp.mp4');
control_vid5 = VideoReader('RHO/296--25--08--(7,15)--42859--RHO-gfp (1).mp4');

% For IL-33, count cells in video over time
num_cells_il33 = cell_counter(il33_vid);
num_cells_il332 = cell_counter(il33_vid2);
num_cells_il333 = cell_counter(il33_vid3);
num_cells_il334 = cell_counter(il33_vid4);
num_cells_il335 = cell_counter(il33_vid5);

% For ORMDL3, count cells in video over time
num_cells_ormdl3 = cell_counter(ormdl3_vid);
num_cells_ormdl32 = cell_counter(ormdl3_vid2);
num_cells_ormdl33 = cell_counter(ormdl3_vid3);
num_cells_ormdl34 = cell_counter(ormdl3_vid4);
num_cells_ormdl35 = cell_counter(ormdl3_vid5);

% For GSDML, count cells in video over time
num_cells_gsdml = cell_counter(gsdml_vid);
num_cells_gsdml2 = cell_counter(gsdml_vid2);
num_cells_gsdml3 = cell_counter(gsdml_vid3);
num_cells_gsdml4 = cell_counter(gsdml_vid4);
num_cells_gsdml5 = cell_counter(gsdml_vid5);

% For TSLP, count cells in video over time
num_cells_tslp = cell_counter(tslp_vid);
num_cells_tslp2 = cell_counter(tslp_vid2);
num_cells_tslp3 = cell_counter(tslp_vid3);
num_cells_tslp4 = cell_counter(tslp_vid4);
num_cells_tslp5 = cell_counter(tslp_vid5);

% For control (RHO gene), count cells in video over time
num_cells_control = cell_counter(control_vid);
num_cells_control2 = cell_counter(control_vid2);
num_cells_control3 = cell_counter(control_vid3);
num_cells_control4 = cell_counter(control_vid4);
num_cells_control5 = cell_counter(control_vid5);

%%
% Fit to unbounded population growth model, assuming population not large
% enough to level off
fit_str = 'c1*exp(a * (x - c2))';
fitmodel = fittype(fit_str);
[fit_ormdl3, metric_ormdl3_bad] = fit((1:length(num_cells_ormdl3))', num_cells_ormdl3', fitmodel);
[fit_gsdml, metric_gsdml_bad] = fit((1:length(num_cells_gsdml))', num_cells_gsdml', fitmodel);
[fit_il33, metric_il33_bad] = fit((1:length(num_cells_il33))', num_cells_il33', fitmodel);
[fit_tslp, metric_tslp_bad] = fit((1:length(num_cells_tslp))', num_cells_tslp', fitmodel);
[fit_control, metric_control_bad] = fit((1:length(num_cells_control))', num_cells_control', fitmodel);

figure
plot(fit_ormdl3, 1:length(num_cells_ormdl3), num_cells_ormdl3)
figure
plot(fit_gsdml, 1:length(num_cells_gsdml), num_cells_gsdml)
figure
plot(fit_il33, 1:length(num_cells_il33), num_cells_il33)
figure
plot(fit_tslp, 1:length(num_cells_tslp), num_cells_tslp)
figure
plot(fit_control, 1:length(num_cells_control), num_cells_control)

% Try fitting to bounded population growth model
fit_str = 'b/(1 + a*exp(-c*x))';
fitmodel = fittype(fit_str);
[fit_ormdl3, metric_ormdl3] = fit((1:length(num_cells_ormdl3))', num_cells_ormdl3', fitmodel, 'Lower', [0 0 0]);
[fit_ormdl32, metric_ormdl32] = fit((1:length(num_cells_ormdl32))', num_cells_ormdl32', fitmodel, 'Lower', [0 0 0]);
[fit_ormdl33, metric_ormdl33] = fit((1:length(num_cells_ormdl33))', num_cells_ormdl33', fitmodel, 'Lower', [0 0 0]);
[fit_ormdl34, metric_ormdl34] = fit((1:length(num_cells_ormdl34))', num_cells_ormdl34', fitmodel, 'Lower', [0 0 0]);
[fit_ormdl35, metric_ormdl35] = fit((1:length(num_cells_ormdl35))', num_cells_ormdl35', fitmodel, 'Lower', [0 0 0]);
[fit_gsdml, metric_gsdml] = fit((1:length(num_cells_gsdml))', num_cells_gsdml', fitmodel, 'Lower', [0 0 0]);
[fit_gsdml2, metric_gsdml2] = fit((1:length(num_cells_gsdml2))', num_cells_gsdml2', fitmodel, 'Lower', [0 0 0]);
[fit_gsdml3, metric_gsdml3] = fit((1:length(num_cells_gsdml3))', num_cells_gsdml3', fitmodel, 'Lower', [0 0 0]);
[fit_gsdml4, metric_gsdml4] = fit((1:length(num_cells_gsdml4))', num_cells_gsdml4', fitmodel, 'Lower', [0 0 0]);
[fit_gsdml5, metric_gsdml5] = fit((1:length(num_cells_gsdml5))', num_cells_gsdml5', fitmodel, 'Lower', [0 0 0]);
[fit_il33, metric_il33] = fit((1:length(num_cells_il33))', num_cells_il33', fitmodel, 'Lower', [0 0 0]);
[fit_il332, metric_il332] = fit((1:length(num_cells_il332))', num_cells_il332', fitmodel, 'Lower', [0 0 0]);
[fit_il333, metric_il333] = fit((1:length(num_cells_il333))', num_cells_il333', fitmodel, 'Lower', [0 0 0]);
[fit_il334, metric_il334] = fit((1:length(num_cells_il334))', num_cells_il334', fitmodel, 'Lower', [0 0 0]);
[fit_il335, metric_il335] = fit((1:length(num_cells_il335))', num_cells_il335', fitmodel, 'Lower', [0 0 0]);
[fit_tslp, metric_tslp] = fit((1:length(num_cells_tslp))', num_cells_tslp', fitmodel, 'Lower', [0 0 0]);
[fit_tslp2, metric_tslp2] = fit((1:length(num_cells_tslp2))', num_cells_tslp2', fitmodel, 'Lower', [0 0 0]);
[fit_tslp3, metric_tslp3] = fit((1:length(num_cells_tslp3))', num_cells_tslp3', fitmodel, 'Lower', [0 0 0]);
[fit_tslp4, metric_tslp4] = fit((1:length(num_cells_tslp4))', num_cells_tslp4', fitmodel, 'Lower', [0 0 0]);
[fit_tslp5, metric_tslp5] = fit((1:length(num_cells_tslp5))', num_cells_tslp5', fitmodel, 'Lower', [0 0 0]);
[fit_control, metric_control] = fit((1:length(num_cells_control))', num_cells_control', fitmodel, 'Lower', [0 0 0]);
[fit_control2, metric_control2] = fit((1:length(num_cells_control2))', num_cells_control2', fitmodel, 'Lower', [0 0 0]);
[fit_control3, metric_control3] = fit((1:length(num_cells_control3))', num_cells_control3', fitmodel, 'Lower', [0 0 0]);
[fit_control4, metric_control4] = fit((1:length(num_cells_control4))', num_cells_control4', fitmodel, 'Lower', [0 0 0]);
[fit_control5, metric_control5] = fit((1:length(num_cells_control5))', num_cells_control5', fitmodel, 'Lower', [0 0 0]);


% Show fits
figure
plot(fit_ormdl32, 1:length(num_cells_ormdl32), num_cells_ormdl32)
figure
plot(fit_gsdml, 1:length(num_cells_gsdml), num_cells_gsdml)
figure
plot(fit_il33, 1:length(num_cells_il33), num_cells_il33)
figure
plot(fit_tslp, 1:length(num_cells_tslp), num_cells_tslp)
figure
plot(fit_control, 1:length(num_cells_control), num_cells_control)

% Aggregate
fit_mat = [fit_ormdl3.c fit_ormdl32.c fit_ormdl33.c fit_ormdl34.c fit_ormdl35.c;
    fit_gsdml.c fit_gsdml2.c fit_gsdml3.c fit_gsdml4.c fit_gsdml5.c;
    fit_il33.c fit_il332.c fit_il333.c fit_il334.c fit_il335.c;
    fit_tslp.c fit_tslp2.c fit_tslp3.c fit_tslp4.c fit_tslp5.c;
    fit_control.c fit_control2.c fit_control3.c fit_control4.c fit_control5.c];

%% Compare rate constants with 1-way ANOVA

data = reshape(fit_mat', 1, 25);
groups = {'ORMDL3', 'ORMDL3', 'ORMDL3', 'ORMDL3', 'ORMDL3', 'GSDML', 'GSDML', 'GSDML', 'GSDML', 'GSDML', 'IL33', 'IL33', 'IL33', 'IL33', 'IL33', 'TSLP', 'TSLP', 'TSLP', 'TSLP', 'TSLP', 'Control', 'Control', 'Control', 'Control', 'Control'};

p = anova1(data, groups)



