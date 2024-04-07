"""
    function read_kml(
        kml_file::String;
        return_table::Bool = false,
        swap_axes::Bool = false,
    )::Table
Reads in a polygon from a KML as a RoamesGeometry Polygon or a Table, currently just finds the first Polygon in the file and reads it.
Also it cannot read polygons with holes at the moment.
"""
function read_kml(
    kml_file::String;
    return_table::Bool = false,
    swap_axes::Bool = false,
)::Union{Table,Polygon}
    s = read(kml_file, String)
    tag = r"<Polygon>"
    m1 = match(tag, s, 1)
    tag = r"<coordinates>"
    m1 = match(tag, s, m1.offset)
    start = m1.offset + length(tag.pattern)
    tag = r"</coordinates>"
    m1 = match(tag, s, m1.offset)
    finish = m1.offset - 1
    temp = split(s[start:finish])

    axis_1 = 2
    axis_2 = 1
    if swap_axes
        axis_1 = 1
        axis_2 = 2
    end
    position = map(temp) do t
        temp2 = split(t, ",")
        if length(temp2) == 2
            SVector(parse(Float64, temp2[axis_1]), parse(Float64, temp2[axis_2]))
        elseif length(temp2) == 3
            SVector(
                parse(Float64, temp2[axis_1]),
                parse(Float64, temp2[axis_2]),
                parse(Float64, temp2[3]),
            )
        end
    end
    position = filter(!isnothing, position)
    return return_table ? Table(position = position) : Polygon(position)
end

"""
    function write_kml(
        edges::polygon::Union{Table,Polygon},
        poly_file::String;
        poly_name::String = "Test",
        color::String = "ffffff",
    )::Nothing
Exports polygon to a `kml` file.
* `edges`: Polygon object or table or edges (positions of the polygon points)
* `poly_file`: Polygon file name
* `poly_name`: Polygon name that goes into the kml file.
* `color`: Color of the polygon

At the moment it only writes single polygon and it does not work with polygons with holes.
"""
function write_kml(
    polygon::Union{Table,Polygon},
    poly_file::String;
    poly_name::String = "Test",
    color::String = "ffffff",
)::Nothing
    edges, has_alt = if polygon isa Table 
        polygon.position, length(polygon.position[1]) == 3
    else
        polygon.exterior.points, length(polygon.exterior.points[1]) == 3
    end
    edges = if has_alt
        map(e -> (lat=e[1], lon=e[2], alt=e[3]), edges)
    else
        map(e -> (lat=e[1], lon=e[2]), edges)
    end
    if is_clockwise(edges)
        edges = reverse(edges)
    end
    open(poly_file, "w") do f
        println(
            f,
            """<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Document>
        <name>$(basename(poly_file))</name>
        <StyleMap id="msn_ylw-pushpin">
                <Pair>
                        <key>normal</key>
                        <styleUrl>#sn_ylw-pushpin</styleUrl>
                </Pair>
                <Pair>
                        <key>highlight</key>
                        <styleUrl>#sh_ylw-pushpin</styleUrl>
                </Pair>
        </StyleMap>
        <Style id="sh_ylw-pushpin">
                <IconStyle>
                        <scale>1.2</scale>
                </IconStyle>
                <LineStyle>
                        <color>ff$color</color>
                </LineStyle>
                <PolyStyle>
                        <color>50$color</color>
                </PolyStyle>
        </Style>
        <Style id="sn_ylw-pushpin">
                <LineStyle>
                        <color>ff$color</color>
                </LineStyle>
                <PolyStyle>
                        <color>50$color</color>
                </PolyStyle>
        </Style>
        <Placemark>
                <name>$(poly_name)</name>
                <styleUrl>#msn_ylw-pushpin</styleUrl>
                <Polygon>
                        <outerBoundaryIs>
                                <LinearRing>
                                        <coordinates>""",
        )
        foreach(edges) do edge_point
            has_alt ? print(f, "$(edge_point.lon),$(edge_point.lat),$(edge_point.alt) ") : print(f, "$(edge_point.lon),$(edge_point.lat) ")
        end
        println(
            f,
            """

                                        </coordinates>
                                </LinearRing>
                        </outerBoundaryIs>
                </Polygon>
        </Placemark>
</Document>
</kml>""",
        )
    end
    return nothing
end