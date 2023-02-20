import React, { useEffect, useRef } from 'react';
import { Container, Row, Col, Form, Button, Card, } from 'react-bootstrap';
import { Jsme } from 'jsme-react';

import ResultsTable from '../Search/ResultsTable';

const url = '/search/similarity';

export default function Similarity() {
    const [smiles, setSmiles] = React.useState(window.location.search.split('smiles=')[1] || '');
    const table = useRef(null);

    useEffect(() => {
        if (window.location.search.split('smiles=')[1]) {
            table.current.getMolecules();
        }
    }, []);

    return (
        <Container fluid className='mt-2'>
            <Row>
                <Col lg={3}>
                    <Card>
                        <Card.Header>Similarity Search</Card.Header>

                        <Card.Body>

                            <p className="mb-3">Search for similar molecules using SMILES.</p>

                            <Form onSubmit={(e) => {
                                e.preventDefault();
                                table.current.getMolecules();
                            }}>
                                <Jsme width="100%" height="300px" smiles={smiles} onChange={(smiles) => setSmiles(smiles)} options="star,newlook" />
                                <Form.Group controlId="smiles" className="mb-3">
                                    <Form.Control type="text" className='mt-1' placeholder="SMILES" value={smiles} onChange={(e) => setSmiles(e.target.value)} />
                                </Form.Group>
                                <Button variant="primary" type="submit">
                                    Search
                                </Button>
                            </Form>
                        </Card.Body>
                    </Card>
                </Col>
                <Col>
                    <ResultsTable url={url} smiles={smiles} ref={table} />
                </Col>
            </Row>
        </Container >
    );
}
