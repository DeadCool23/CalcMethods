is_find_max = false;
is_debug = true;

mtx = [
    7 4 3 8 2;
    4 5 1 6 3;
    8 4 5 7 2;
    1 2 4 7 2;
    3 9 9 2 5
];

venger_method(mtx, is_find_max, is_debug)

function venger_method(srcC, is_find_max, is_debug)
    fprintf('Матрица стоимостей:\n');
    C = srcC;
    disp(C);

    if is_find_max
        C = reduction_to_minimization(C, is_debug);
    end

    [C, SIZ] = prepare_step(C, is_debug);

    X_opt = main_step(C, SIZ, is_debug);
    fprintf("X_opt = \n");
    disp(X_opt);

    f_opt = get_f_opt(srcC, X_opt);
    fprintf("f_opt = %d\n", f_opt)
end

function C = reduction_to_minimization(C, is_debug)
    max_el = max(max(C));
    C = max_el + -1 * C;

    if is_debug
        fprintf("Максимальный элемент: %d\n", max_el);
        fprintf("Сведенная к задаче минимизации матрица стоимостей:\n");
        disp(C);
    end
end

function [C, SIZ] = prepare_step(C, is_debug)
    if is_debug
        fprintf("Подготовительный этап:\n\n");
    end

    C = sub_min_cols(C, is_debug);
    C = sub_min_rows(C, is_debug);
    SIZ = init_SIZ(C, is_debug);
end

function X_opt = main_step(C, SIZ, is_debug)
    if is_debug
        fprintf("\nОсновной этап:\n\n");
    end

    k = get_SIZ_k(SIZ);
    [n, ~] = size (C);

    iters = n - k;
    if is_debug && iters ~= 0
        fprintf("k = %d < n = %d. Требуется %d итераций\n\n", k, n, iters)
    elseif is_debug
        fprintf("k = %d == n = %d. Итераций не требуется\n\n", k, n)
    end
    
    iter_num = 1;
    while iters ~= 0
        if is_debug
            fprintf("Итерация №%d\n\n", iter_num);
        end

        [C, SIZ] = improvement_SIZ(C, SIZ, is_debug);

        iter_num = iter_num + 1;
        iters = iters - 1;
    end

    X_opt = SIZ;
end

function mtx = sub_min_cols(mtx, is_debug)
    min_els = min(mtx);
    for i = 1 : length(min_els)
        mtx(:, i) = mtx(:, i) - min_els(i);
    end
    
    if is_debug
        fprintf("Минимальные элементы столбцов:\n");
        disp(min_els);
        fprintf("Матрица стоимостей после вычитания минимальных элементов из столбцов:\n");
        disp(mtx);
    end
end

function mtx = sub_min_rows(mtx, is_debug)
    min_els = min(mtx, [], 2);
    for i = 1 : length(min_els)
        mtx(i, :) = mtx(i, :) - min_els(i);
    end

    if is_debug
        fprintf("Минимальные элементы строк:\n");
        disp(min_els);
        fprintf("Матрица стоимостей после вычитания минимальных элементов из строк:\n");
        disp(mtx);
    end
end

function SIZ = init_SIZ(C, is_debug)
    [rows, cols] = size(C);
    SIZ = zeros(rows, cols);
    
    for i = 1 : rows
        for j = 1 : cols
            if C(i,j) == 0 && ~any(SIZ(:, j) == 1)
                SIZ(i, j) = 1;
                break
            end
        end
    end

    if is_debug
        fprintf("Начальная СНН:\n");
        print_SIZ(C, SIZ);
        fprintf("Количество нулей СНН: k = %d\n", get_SIZ_k(SIZ));
    end
end

function [C, SIZ] = improvement_SIZ(C, SIZ, is_debug)
    [rows, cols] = size(SIZ);
    marked_zeroes = zeros(rows, cols);
    selected_cols = sum(SIZ);
    selected_rows = zeros(rows);
    selected = update_selected(rows, cols, selected_rows, selected_cols);

    if is_debug
        fprintf("СНН с начально выделенными столбцами:\n");
        print_marked_SIZ(C, SIZ, marked_zeroes, selected_rows, selected_cols)
    end

    marked_zero_pos = [-1 -1];
    is_L_chain_builded = false;
    while ~is_L_chain_builded
        marked_zero_pos = get_zero_mark_pos(C, selected);
        
        if marked_zero_pos(1) == -1
            C = add_zeros_to_unselected(C, selected, selected_rows, selected_cols);

            if is_debug
                fprintf('\nСреди невыделенных элементов нет 0, преобразуем матрицу:\n');
                print_marked_SIZ(C, SIZ, marked_zeroes, selected_rows, selected_cols);
            end

            marked_zero_pos = get_zero_mark_pos(C, selected);
        end
        
        marked_zeroes = mark_zero(marked_zeroes, marked_zero_pos);
        if is_debug
            fprintf("\nОтметим 0' найденный невыделенный 0:\n");
            print_marked_SIZ(C, SIZ, marked_zeroes, selected_rows, selected_cols);
        end

        stared_zero_pos = get_stared_zero_pos_on_row_with_marked(SIZ, marked_zero_pos);
        if stared_zero_pos(1) == -1
            is_L_chain_builded = true;
        else
            selected_rows(stared_zero_pos(1)) = 1;
            selected_cols(stared_zero_pos(2)) = 0;
            selected = update_selected(rows, cols, selected_rows, selected_cols);

            if is_debug
                fprintf("\nПеренесем выделение со столбца на строку, т.к. на строке с 0' есть 0*:\n");
                print_marked_SIZ(C, SIZ, marked_zeroes, selected_rows, selected_cols);
            end
        end
    end

    SIZ = build_L_chain(SIZ, marked_zeroes, marked_zero_pos, is_debug);
end

function SIZ = build_L_chain(SIZ, marked_zeros, marked_zero_pos, is_debug)
    if is_debug
        fprintf("\nПостроена L-цепочка: ")
    end
    
    i = marked_zero_pos(1);
    j = marked_zero_pos(2);

    is_L_chain_ended = false;
    while ~is_L_chain_ended
        marked_zeros(i, j) = 0;

        stared_zero_row_pos = find(SIZ(:, j) == 1, 1);
        SIZ(i, j) = 1;

        if is_debug
            fprintf("0' (%d, %d) ", i, j);
        end
        
        if isempty(stared_zero_row_pos)
            is_L_chain_ended = true;
        else
            SIZ(stared_zero_row_pos, j) = 0;

            if is_debug
                fprintf("-> 0* (%d, %d) -> ", stared_zero_row_pos, j);
            end

            marked_zero_col_pos = find(marked_zeros(stared_zero_row_pos, :) == 1, 1);
            i = stared_zero_row_pos;
            j = marked_zero_col_pos;
        end
    end

    if is_debug
        fprintf("\n\n");
    end
end

function marked_zeroes = mark_zero(marked_zeroes, marked_zero_pos)
    marked_zeroes(marked_zero_pos(1), marked_zero_pos(2)) = 1;
end

function C = add_zeros_to_unselected(C, selected, selected_rows, selected_cols)
    h = find_h(C, selected);
    [rows, cols] = size(C);

    for i = 1 : cols
        if selected_cols(i) == 0
            C(:, i) = C(:, i) - h;
        end 
    end
    for i = 1 : rows
        if selected_rows(i) == 1
            C(i, :) = C(i, :) + h;
        end 
    end
end

function h = find_h(C, selected)
    h = inf;
    [rows, cols] = size(C);
    for i = 1 : rows
        for j = 1 : cols
            if selected(j, i) == 0 && C(j, i) < h
                h = C(j, i);
            end
        end 
    end
end

function marked_zero_pos = get_zero_mark_pos(C, selected)
    marked_zero_pos = [-1 -1];
    [rows, cols] = size(C);
    for j = 1 : cols
        for i = 1 : rows
            if selected(i, j) == 0 && C(i, j) == 0
                marked_zero_pos(1) = i;
                marked_zero_pos(2) = j;
                return
            end
        end
    end
end

function stared_zero_pos = get_stared_zero_pos_on_row_with_marked(SIZ, marked_zero_pos)
    i = marked_zero_pos(1);
    [~, cols] = size(SIZ);

    stared_zero_pos = [-1 -1];
    for j = 1 : cols
       if SIZ(i, j) == 1
           stared_zero_pos(1) = i;
           stared_zero_pos(2) = j;
           break
       end 
    end
end

function selected = update_selected(rows, cols, selected_rows, selected_cols)
    selected = zeros(rows, cols);
    selected(:, selected_cols == 1) = 1;
    selected(selected_rows == 1, :) = 1;
end

function k = get_SIZ_k(SIZ)
    k = sum(sum(SIZ));
end

function print_SIZ(C, SIZ)
    [rows, cols] = size(C);
    for i = 1 : rows
        for j = 1 : cols
            if SIZ(i, j) == 1
                fprintf("\t%d*", C(i,j));
            else
                fprintf("\t%d", C(i,j));
            end
        end
        fprintf("\n");
    end 
end

function print_marked_SIZ(C, SIZ, marked_zeros, selected_rows, selected_cols)
    [rows, cols] = size(C);
    for i = 1 : rows
        for j = 1 : cols
            if SIZ(i, j) == 1
                fprintf("\t%d*", C(i,j));
            elseif marked_zeros(i, j) == 1
                fprintf("\t%d'", C(i,j));
            else
                fprintf("\t%d", C(i,j));
            end
        end
        if selected_rows(i) == 1
            fprintf("\t+");
        end
        fprintf("\n");
    end 

    for i = 1 : cols
        if selected_cols(i) == 1
            fprintf("\t+");
        else
            fprintf("\t ");
        end
    end
    fprintf("\n");
end

function f_opt = get_f_opt(C, X_opt)
    f_opt = sum(sum(C .* X_opt));
end