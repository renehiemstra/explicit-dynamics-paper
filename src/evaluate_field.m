function f = evaluate_field(field, granularity)

    % initialize
    u.p = field.basis{1}.p; u.kts = field.basis{1}.kts;
    v.p = field.basis{2}.p; v.kts = field.basis{2}.kts;

    % compute basis functions
    b.u = field.basis{1}.eval.trialfuns(linspace(u.kts(1), u.kts(end), granularity(1))');
    b.v = field.basis{2}.eval.trialfuns(linspace(v.kts(1), v.kts(end), granularity(2))');

    % compute coordinates
    f = b.u{1} * field.coeffs * b.v{1}';
end