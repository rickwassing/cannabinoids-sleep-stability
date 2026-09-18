function str = stringify_pvalue(p)
if p < 1e-4
    str = sprintf('<10^%.0f', round(log10(p) + 1));
elseif p < 0.001
    str = sprintf('%.5f', p);
elseif p < 0.01
    str = sprintf('%.4f', p);
else
    str = sprintf('%.3f', p);
end
end
