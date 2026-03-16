% Create unique identifier
function id = getuuid()
id = char(matlab.lang.internal.uuid());
id = id(1:8);
end