load outmat/specparam.mat

[pgrid(:,1:4) vdgrid cd1grid cd2grid ce1grid ce2grid]

[pgrid(:,1:4) cd2grid ce2grid b*cd1grid.*ce2grid-ce1grid cd1grid.*(b-cd2grid)]