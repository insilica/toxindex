CREATE TABLE chemical_properties (
    chemical_property_id SERIAL PRIMARY KEY,
    chemical_id INTEGER REFERENCES chemicals(chemical_id),
    property_id INTEGER REFERENCES properties(property_id),
    value TEXT
);