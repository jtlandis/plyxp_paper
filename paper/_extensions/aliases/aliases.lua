-- I do not believe the panache LSP will
-- auto complete these short code :(
return {
    ['plyxp'] = function(args, kwargs, meta, raw_args, context)
        -- see https://quarto.org/docs/extensions/shortcodes.html
        -- for documentation on shortcode development
        return pandoc.Emph(pandoc.Str("plyxp"))
    end,
    ['dplyr'] = function(args, kwargs, meta, raw_args, context)
        return pandoc.Emph(pandoc.Str("dplyr"))
    end,
    ['purrr'] = function(args, kwargs, meta, raw_args, context)
        return pandoc.Emph(pandoc.Str("purrr"))
    end,
    ['pandas'] = function(args, kwargs, meta, raw_args, context)
        return pandoc.Emph(pandoc.Str("pandas"))
    end,
    ['anndata'] = function(args, kwargs, meta, raw_args, context)
        return pandoc.Emph(pandoc.Str("AnnData"))
    end,
    ['rlang'] = function(args, kwargs, meta, raw_args, context)
        return pandoc.Emph(pandoc.Str("rlang"))
    end
}
