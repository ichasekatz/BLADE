-- Give every table full page width, with column widths proportional to the longest cell text (clamped), so that
-- short pipe tables are not squeezed into narrow auto-width columns in DOCX.
local function cell_len(cell)
  return #pandoc.utils.stringify(cell.contents)
end

function Table(tbl)
  local n = #tbl.colspecs
  local maxlen = {}
  for i = 1, n do maxlen[i] = 7 end
  local function scan(rows, header)
    for _, row in ipairs(rows) do
      for i, cell in ipairs(row.cells) do
        if i <= n then
          if not header then maxlen[i] = math.max(maxlen[i], math.min(cell_len(cell), 40)) end
          for w in pandoc.utils.stringify(cell.contents):gmatch("%S+") do maxlen[i] = math.max(maxlen[i], math.min(#w + 3, 22)) end
        end
      end
    end
  end
  scan(tbl.head.rows, true)
  for _, body in ipairs(tbl.bodies) do scan(body.body, false) end
  local total = 0
  for i = 1, n do total = total + maxlen[i] end
  for i = 1, n do
    tbl.colspecs[i] = {tbl.colspecs[i][1], maxlen[i] / total}
  end
  return tbl
end
