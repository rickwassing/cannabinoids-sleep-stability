function str = parsepvalue(p)

if p < 0.0001
    str = sprintf('<10^%.0f', round(log10(p)+1));
elseif p < 0.001
    str = sprintf('%.5f', p);
elseif p < 0.01
    str = sprintf('%.4f', p);
elseif p < 0.05
    str = sprintf('%.3f', p);
else
    str = sprintf('%.3f', p);
end

end