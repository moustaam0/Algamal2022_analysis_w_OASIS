Project_path = '/Users/moustafaalgamal/Documents/OASIS_matlab-master';
name = 'P54F4_A';

smin = -2.5; % minimum event size relevant to noise standard deviation(-2 = two times SD)
lambda = 0; % regularaization parameter for the model (set to 0 and increase smin, if unkwon)
thr_inactive = 0.001; % theshould for inactive cells

analyzeCa1 (Project_path, name);
corrfix_final3(Project_path, name);
