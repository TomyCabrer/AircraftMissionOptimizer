function MissionPLot(cp, mission, opts, out)
%MISSIONPLOT  Event-based mission profile plot (altitude + weight fraction).
%
% Purpose
%   Visualise the design mission as two subplots:
%     1) Altitude at key mission nodes (0..9)
%     2) Cumulative weight fraction W/W0 across the same nodes
%
% Inputs
%   cp      : constants struct (contains fixed segment fractions: W1warm, W2climb, W4dec, W10land)
%   mission : mission struct (used here for loiter altitude, etc.)
%   opts    : solver options (not used directly, kept for a consistent API)
%   out     : solver output struct (from MissionProfileSolver)
%
% Notes
%   - Plot is "event-based": the x-axis is node index, not time or distance.
%   - Points are labelled with readable markers + numeric annotations.
%   - Replace placeholder fractions WF(5) and WF(8) if you later model those segments.

% -------------------- Unpack mission (only what is needed) ----------------
[~, ~, ~, ~, ~, ~, ~, ~, h_loiter, ~, ~, ~, ~, ~, ~, ~, ~] = unpackMissionOpts(mission, opts);

% -------------------- Unpack solver output (only what is needed) ----------
[ ...
    ~, ~, ...
    hcruise_start, hcruise_end, hcruise2_start, hcruise2_end, ...
    ~, ~, ...
    ~, ~, ...
    ~, ~, ...
    ~, ~, ...
    ~, ~, ...
    W1cruise, W2cruise, Wloiter, ...
    ~, W7decloit, W5climb, ...
    ~, ~, ~, ~, ...
    ~, ~, ~, ...
    ~, ~, ...
    ~, ~, ...
    ~, ~, ...
    ~, ~, ~, ~ ] = unpackOut(out);

% -------------------- Mission definition (nodes & labels) -----------------
seg = ["0-1 Taxi/TO", "1-2 Climb+accel", "2-3 Cruise (main)", "3-4 Descent", ...
       "4-5 Missed+climb", "5-6 Cruise (alt)", "6-7 Loiter", "7-8 Descent", "8-9 Land+Taxi"];
x = 0:9;

% Altitude at mission nodes (0..9)
h_node = [ ...
    0;              % 0 start
    0;              % 1 after taxi/TO
    hcruise_start;  % 2 end climb to main cruise
    hcruise_end;    % 3 end main cruise
    0;              % 4 after descent (field)
    hcruise2_start; % 5 end missed+climb to alternate cruise
    hcruise2_end;   % 6 end alternate cruise
    h_loiter;       % 7 loiter altitude (from mission input)
    0;              % 8 after descent
    0               % 9 after landing+taxi
];

% -------------------- Weight fractions per segment ------------------------
% WF(k) = W_after/W_before for segment (k-1 -> k), k=1..9
WF = ones(1,9);
WF(1) = cp.W1warm;     % 0-1 Taxi/TO + warmup
WF(2) = cp.W2climb;    % 1-2 Climb to cruise 1
WF(3) = W1cruise;      % 2-3 Cruise 1
WF(4) = cp.W4dec;      % 3-4 Descent/decel
WF(5) = W5climb;           % 4-5 Missed+climb (placeholder)
WF(6) = W2cruise;      % 5-6 Cruise 2
WF(7) = Wloiter;       % 6-7 Loiter
WF(8) =  W7decloit;           % 7-8 Descent (placeholder)
WF(9) = cp.W10land;    % 8-9 Landing/taxi

% Cumulative weight fraction at each node (W/W0)
Wfrac_node = ones(size(x));
for k = 2:numel(x)
    Wfrac_node(k) = Wfrac_node(k-1) * WF(k-1);
end

% -------------------- Plot settings for readability -----------------------
ms = 7;      % marker size
lw = 1.6;    % line width
fs = 10;     % base font size
fs_anno = 9; % annotation font size

figure("Name","Mission Profile","Color","w");

% ===================== Subplot 1: Altitude =====================
subplot(2,1,1);
plot(x, h_node, "-o", "LineWidth", lw, "MarkerSize", ms);
grid on;
ylabel("Altitude [m]");
title("Mission profile (event-based)");
xticks(x);
xticklabels(string(x));
xlim([0 9]);
set(gca,"FontSize",fs);

% Alternate label placement (up / down) + background box
dy_h = 0.04 * max(h_node);     % vertical spacing
for i = 1:numel(x)
    if mod(i,2)==0
        yoff = +dy_h;
    else
        yoff = -dy_h;
    end
    text(x(i), h_node(i) + yoff, sprintf("%.0f m", h_node(i)), ...
        "FontSize", fs_anno, ...
        "HorizontalAlignment","center", ...
        "BackgroundColor","w", ...
        "EdgeColor",[0.6 0.6 0.6], ...
        "Margin",2);
end

% ===================== Subplot 2: Weight fraction =====================
subplot(2,1,2);
plot(x, Wfrac_node, "-o", "LineWidth", lw, "MarkerSize", ms);
grid on;
ylabel("Weight fraction W/W_0 [-]");
xlabel("Mission node index");
xticks(x);
xticklabels(string(x));
xlim([0 9]);
set(gca,"FontSize",fs);

% Label each point with weight fraction value
dy_w = 0.02 * (max(Wfrac_node) - min(Wfrac_node) + eps);
for i = 1:numel(x)
    text(x(i) + 0.06, Wfrac_node(i) + dy_w, sprintf("%.4f", Wfrac_node(i)), ...
        "FontSize", fs_anno, "Interpreter","none");
end


for k = 1:9
    xm = k - 0.5;
    ym = 0.5 * (Wfrac_node(k) + Wfrac_node(k+1));
    text(xm, ym, sprintf("WF=%.4f", WF(k)), ...
        "FontSize", fs_anno, ...
        "HorizontalAlignment","center", ...
        "BackgroundColor","w", ...
        "EdgeColor",[0.7 0.7 0.7], ...
        "Margin",2);
end

% Optional: add a compact legend-style text box with segment names
% (kept short to avoid clutter; comment out if undesired)
annotation("textbox", [0.13 0.92 0.85 0.06], "String", strjoin(seg, " | "), ...
    "EdgeColor", "none", "FontSize", 8, "Interpreter","none");

end