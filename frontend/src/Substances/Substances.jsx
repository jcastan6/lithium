import React, { useEffect, useRef } from 'react';
import { Container, Row, Col, Form, Button, Card, } from 'react-bootstrap';
import { Jsme } from 'jsme-react';

import SubstanceTable from '../Search/SubstanceTable';

const url = '/substance';


export default function Substances() {
    const [smiles, setSmiles] = React.useState(window.location.search.split('smiles=')[1] || '');
    const table = useRef(null);

    useEffect(() => {
        table.current.getMolecules();
    }, []);

    return (
        <Container fluid className='mt-2'>
            <Row>
                <Col>
                    <SubstanceTable url={url} ref={table} />
                </Col>
            </Row>
        </Container >
    );
}
