function ang = circ_rad2deg360(alpha)
% Converts radians in [-pi pi] to degrees in [0 360)
% 0° (Peak) — 90° (Descending) — 180° (Trough) — 270° (Ascending) — 360° (Peak)

ang = mod(alpha, 2*pi) * 180/pi;
end