#!/usr/bin/env julia

"""
    md_to_jl(mdfile, jlfile)

Convert a Markdown file to a Julia file:

- Fenced code blocks (```...``` possibly with an environment) → lines inside them become verbatim code.
- Lines outside code blocks → commented with `#`.
- Fence lines become a special marker line `##!CODEBOUNDARY @example SomeEnv`.
"""
function md_to_jl(mdfile::String, jlfile::String)
    inside_code_block = false
    output_lines = String[]
    
    # Regex to match a line that starts with triple backticks (optionally with environment text)
    # Examples matched: ```     ```julia     ```@example test
    rx_fence = r"^```+(\s*(.*))?$"

    for line_raw in eachline(mdfile)
        line = replace(chomp(line_raw), '\ufeff' => "")  # strip BOM, trailing CR

        # Check if line is a fence start/end
        m = match(rx_fence, line)
        if m !== nothing
            # Toggle code-block mode
            inside_code_block = !inside_code_block
            # Extract environment text, e.g. "julia", "@example test", or empty
            env = m.captures[2] === nothing ? "" : strip(m.captures[2])
            if env == ""
              push!(output_lines, "##!CODEBOUNDARY")
              push!(output_lines, "")
            else
              push!(output_lines, "")
              push!(output_lines, "##!CODEBOUNDARY $env")
            end
            continue
        end

        if inside_code_block
            # Inside fenced code → write verbatim
            push!(output_lines, line)
        else
            # Outside fenced code → comment it
            if isempty(strip(line))
                push!(output_lines, "#")
            else
                push!(output_lines, "# " * line)
            end
        end
    end

    # Write to output file
    open(jlfile, "w") do f
        for l in output_lines
            println(f, l)
        end
    end
end

"""
    jl_to_md(jlfile, mdfile)

Convert a Julia file (produced by `md_to_jl`) back to Markdown:

- Lines like `##!CODEBOUNDARY something` toggle opening/closing a fenced code block with
  ```````something````` (if `something` is nonempty).
- Lines starting with `#` outside code blocks become normal text (the `#` is stripped).
- Other lines become code lines (inside code blocks).
- Empty lines inside code blocks remain code; empty lines outside remain text (no extra fences).
"""
function jl_to_md(jlfile::String, mdfile::String)
    output_lines = String[]
    code_block_open = false

    # Regex to detect boundary marker lines, e.g. "##!CODEBOUNDARY @example"
    rx_marker = r"^##!CODEBOUNDARY\s*(.*)$"

    # Helper to close any open code block in the output
    function close_block!()
        if code_block_open
            push!(output_lines, "```")
            code_block_open = false
        end
    end

    open(jlfile, "r") do f
        for line_raw in eachline(f)
            line = replace(chomp(line_raw), '\ufeff' => "")  # strip BOM, trailing CR

            # Check for code-boundary marker
            m = match(rx_marker, line)
            if m !== nothing
                # We have a boundary line like "##!CODEBOUNDARY someenv"
                env = strip(m.captures[1])  # might be empty or "julia", "@example XY", etc.

                if code_block_open
                    # If already in a code block, close it
                    push!(output_lines, "```")
                    code_block_open = false
                else
                    # Otherwise, open a code block with the environment
                    if isempty(env)
                        push!(output_lines, "```")
                    else
                        push!(output_lines, "```" * env)
                    end
                    code_block_open = true
                end
                continue
            end

            # Not a boundary marker:
            if code_block_open
                # We are inside a code block → all lines go verbatim as code,
                # even if they're empty. So just push it.
                push!(output_lines, line)

            else
                # We're outside a code block → treat lines as text unless we have
                # a reason to open a code block.

                # If line starts with '#' → strip '#' and treat as text
                if startswith(strip(line), "#")
                    text_line = strip(strip(line, ['#']))
                    push!(output_lines, text_line)

                else
                    # If it's an empty line outside code blocks, treat it as text
                    # so we don't spuriously open a code block just for an empty line.
                    if isempty(strip(line))
                        push!(output_lines, "")
                    else
                        # Non-empty, non-`#` line → code outside a block → open a block
                        push!(output_lines, "```")
                        code_block_open = true
                        push!(output_lines, line)
                    end
                end
            end
        end
    end

    # If we're still in a code block at EOF, close it
    close_block!()

    # Write final .md
    open(mdfile, "w") do f
        for l in output_lines
            println(f, l)
        end
    end
end

# ------------------------------------------------------------------------
# Command-line entry point
# ------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 3
        println("Usage:")
        println("  julia md_jl_converter.jl md_to_jl <input.md> <output.jl>")
        println("  julia md_jl_converter.jl jl_to_md <input.jl> <output.md>")
        exit(1)
    end

    mode = ARGS[1]
    fin  = ARGS[2]
    fout = ARGS[3]

    if mode == "md_to_jl"
        md_to_jl(fin, fout)
    elseif mode == "jl_to_md"
        jl_to_md(fin, fout)
    else
        println("Unknown mode: $mode")
        exit(1)
    end
end
