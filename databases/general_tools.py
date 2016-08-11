import sqlalchemy


def get_or_create(model, session, **kwargs):
    """Uses kwargs to get instance of model. If can't get that instance from session, add it to session.

    Notes
    -----
    This is an analogue to the Django function get_or_create, written for sqlalchemy.
    Implements sqlalchemy query function one_or_none(). See:
        http://docs.sqlalchemy.org/en/latest/orm/query.html#sqlalchemy.orm.query.Query.one_or_none
    Code adapted from a response to a question on Stack Overflow:
        http://stackoverflow.com/questions/2546207/does-sqlalchemy-have-an-equivalent-of-djangos-get-or-create
    Posted by stackoverflow.com user Wolph:
        http://stackoverflow.com/users/54017/wolph

    Parameters
    ----------
    session : session
        An sqlalchemy session.
    model : Base
        An sqlalchemy table.
    kwargs : dict or specified keyword=value pairs
        key/value pairs used to instantiate the model.

    Returns
    -------
    2-tuple:
        First element of tuple is instance of model instantiated by **kwargs
        Second element of tuple is bool:
            True if instance has been added to session.
            False is instance was already present.

    Raises
    ------
    sqlalchemy.orm.exc.MultipleResultsFound
        If query selects multiple objects.
    """
    instance = session.query(model).filter_by(**kwargs).one_or_none()
    if instance:
        return instance, False
    else:
        params = dict((k, v) for k, v in kwargs.items() if not isinstance(v, sqlalchemy.sql.expression.ClauseElement))
        instance = model(**params)
        session.add(instance)
        # Need this step to retrieve automatically populated columns, i.e., the assigned id
        instance = session.query(model).filter_by(**kwargs).one()
        return instance, True


def populate_model(model, data_dicts, session):
    """Create instances of model with data from data_dicts, and add them to session.

    Runs get_or_create on the model, using the data in each data_dict.

    Parameters
    ----------
    model : Base
        An sqlalchemy table.
    data_dicts : list of dict
        Each dictionary represents the attributes to instantiate the model.
        Keys in the dictionaries are the names of model fields.
        Vals are the data to be inserted into the row corresponding to the field.
    session : An sqlalchemy session.

    Returns
    -------
    created_objs : list, or None
        List of model objects created, if any.
    """
    created_objs = []
    for data_dict in data_dicts:
        obj, created = get_or_create(model=model, session=session, **data_dict)
        if created:
            created_objs.append(obj)
    return created_objs


__author__ = 'Kieran L. Hudson, Jack W. Heal'
